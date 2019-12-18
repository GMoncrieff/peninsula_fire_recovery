///////////////////////////////////////
//Fire recovery anomalies in fynbos
//Jasper Slingsby, Glenn Moncrieff and Adam Wilson
////////////////////////////////////////
var batch = require('users/fitoprincipe/geetools:batch');

//how long to monitor
var length = -6;
var aoi = 
    ee.Geometry.Polygon(
        [[[18.158159615255272, -33.861166832966596],
          [18.158159615255272, -34.41584793013852],
          [18.740435005880272, -34.41584793013852],
          [18.740435005880272, -33.861166832966596]]], null, false);
var aoiras = ee.FeatureCollection(ee.Feature(aoi,{flag:1}))
  .reduceToImage({
    properties: ['flag'],
    reducer: ee.Reducer.first()
  })

//cloud/quality masking
var getQABits = function(image, start, end, newName) {
    // Compute the bits we need to extract.
    var pattern = 0;
    for (var i = start; i <= end; i++) {
       pattern += Math.pow(2, i);
    }
    return image.select([0], [newName])
                  .bitwiseAnd(pattern)
                  .rightShift(start);
};

var qualitymask = function(image) {
  // Select the QA band
  var QA = image.select('DetailedQA');
  
  // Get the cloud_state bits and find cloudy areas.
  var cloud = getQABits(QA, 0, 1, 'MODLAND')
                      .expression("b(0) < 2");
                      
  return(image.updateMask(cloud))
  
};
  
// phi parameter
var phi = -1.56;
// sigma parameter
var sigma = 0.053;
//get today date
var today = ee.Date(new Date()).millis();
//apply algorthym over this period (x months)
var start = ee.Date(new Date()).advance(length, 'month');


//get recent ndvi data
var mod_index = ee.ImageCollection('MODIS/006/MYD13Q1')
.merge(ee.ImageCollection('MODIS/006/MOD13Q1'))
.filter(ee.Filter.date(start, today))
.filterBounds(aoi)
.map(function(image) {
  return image.updateMask(aoiras)
})
.map(qualitymask)
.select('NDVI')
.map(function(image) {
  return image.addBands(image.select('NDVI').divide(10000));
})
.select(['NDVI_1'],['NDVI']);


/////////////////
//time since fire
//////////////////

//function to get time since fire (nested to pass multiple vars)
//https://gis.stackexchange.com/questions/302760/imagecollection-map-with-multiple-input-function

var getTSB_o = function(idate) {
  var getTSB_i  = function(image){
    var year = ee.Date(image.get('system:time_start')).get('year');
    var st_date = ee.Date(ee.String(year).cat('-01-01')).millis();
    return image.multiply(86400000).add(st_date).subtract(idate).abs().divide(86400000);
  }
return getTSB_i;
}

//function to add time since fire as a band
function addTSB(image){
  var idate = ee.Date(image.get('system:time_start')).millis();
  var tsb_col = ee.ImageCollection("MODIS/006/MCD64A1")
  .filter(ee.Filter.date('2000-01-01',idate))
  .filterBounds(aoi)
  .map(function(image) {
  return image.updateMask(aoiras)
  })
  .select(['BurnDate'])
  .map(getTSB_o(idate));
  var tsb_min = tsb_col.min();
  return image.addBands(tsb_min);
}

//get the month of each fire
function firemonth(image){
  //date in mill since 1970
  var bdate = ee.Image(ee.Date(image.get('system:time_start')).millis()).rename(['cDate']);
  //age in mill
  var bdate = bdate.addBands(image.select(['BurnDate']).multiply(86400000).rename(['BurnDateSec']))
  //calc burn date
  var bdate = bdate.addBands(bdate.select(['cDate']).subtract(bdate.select(['BurnDateSec'])).rename(['cnDate']))
  //calc burn month
  var bmonth = bdate.select('cnDate').mod(ee.Image(31557600000)).divide(2629800000).ceil().rename(['Month']);
  
  return image.addBands(bmonth);
}

//calculate fire age and add to collection
var mod_data = mod_index
.map(addTSB)
.map(firemonth);


/////////////////////
//add noise to pars
////////////////////

//number of samples from normal dist
var sample = 100;

//read in data par a
var amean = ee.Image("users/glennwithtwons/peninsulapars/alphaM");
var asd = ee.Image("users/glennwithtwons/peninsulapars/alphaSD");

//read in data par a
var bamean = ee.Image("users/glennwithtwons/peninsulapars/AM");
var basd = ee.Image("users/glennwithtwons/peninsulapars/ASD");

//read in data par gamma
var gmean =  ee.Image("users/glennwithtwons/peninsulapars/gammaM");
var gsd = ee.Image("users/glennwithtwons/peninsulapars/gammaSD");

//read in data par lambda
var lmean = ee.Image("users/glennwithtwons/peninsulapars/lambdaM");
var lsd = ee.Image("users/glennwithtwons/peninsulapars/lambdaSD");

//var lmean = ee.Image.constant(10);
//var lsd = ee.Image.constant(5);

//create normal random var
var rnormal = function(seed,step,mean,sd) {
  //create unif randoms
  var v = ee.Image.random(seed);
  var u = ee.Image.random(seed.add(step));
  var r = u.log().multiply(-2).sqrt().multiply(v.multiply(2*Math.PI).cos());
  return r.multiply(sd).add(mean);
};

//create normal image samples (nested to pass multiple vars)
//https://gis.stackexchange.com/questions/302760/imagecollection-map-with-multiple-input-function

//outer function
var norm_rep_o = function(image){
//inner function
  var norm_rep_i  = function(s) {
    var ns = ee.Number(s);
    
    //calc normal
    var anorm = rnormal(ns,sample,amean,asd).select(['random'],['anorm']);
    var gnorm = rnormal(ns.add(sample),sample,gmean,gsd).select(['random'],['gnorm']);
    var lnorm = rnormal(ns.add(sample).add(sample),sample,lmean,lsd).select(['random'],['lnorm']);
    var banorm = rnormal(ns.add(sample).add(sample).add(sample),sample,bamean,basd).select(['random'],['banorm']);
    
    //fit curve
    var acalc = image.addBands(anorm).addBands(gnorm).addBands(lnorm).addBands(banorm);
    acalc = acalc.addBands(acalc.select(['BurnDate'],['expon']).divide(365.25).divide(acalc.select('lnorm')).multiply(-1).exp());
    acalc = acalc.addBands(acalc.select(['Month'],['sinpar']).subtract(1).multiply(3.141593/6).add(phi)
    .add(acalc.select(['BurnDate']).divide(365.25).multiply(6.283185)).sin().multiply(acalc.select('banorm')));
    
    
    //calc expected ndvi
    var obs_mean = acalc.expression(
      'ALPHA+GAMMA-GAMMA*EXPO+SINPAR',{
        'BD': acalc.select('BurnDate'),
        'ALPHA': acalc.select('anorm'),
        'GAMMA': acalc.select('gnorm'),
        'EXPO': acalc.select('expon'),
        'SINPAR': acalc.select('sinpar')
    })
    
    var obs_norm = rnormal(ns.add(sample).add(sample).add(sample).add(sample),sample,obs_mean,sigma).select(['random'],['banorm']);
    return(obs_norm.rename(['anorm']).clamp(0,1));
  };
  return norm_rep_i;
};

//apply to each each in ndvi/age collection
var summary_norm = function(image){
  //create random seed seq
  var tseed = ee.Date(new Date()).millis();
  var seq = ee.List.sequence(tseed,tseed.add(sample),1);
  //get list of image
  var imlist = seq.map(norm_rep_o(image));
  //create an image coll of expected_ndvi samples and calc percentiles
  var percent = ee.ImageCollection(imlist).reduce(ee.Reducer.percentile([5,95]));
 
  return image.addBands(percent);
};

//calculate percentiles
var ndvi_percent = mod_data.map(summary_norm);
//print(ndvi_percent)

//have percentiles been exceeded
var exceed_a = ndvi_percent.map(function(image) {
  return image.expression(
      '(NDVI > B95)',{ 
        'NDVI': image.select('NDVI'),
        'B95': image.select('anorm_p95')
    }).select([0],['exceed_above'])
    .set('system:time_start',image.get('system:time_start'));
});

var exceed_b = ndvi_percent.map(function(image) {
  return image.expression(
      '(NDVI < B5)',{ 
        'NDVI': image.select('NDVI'),
        'B5': image.select('anorm_p5')
    }).select([0],['exceed_below'])
    .set('system:time_start',image.get('system:time_start'));
});

////////////////////////////
//summarise exceedances
////////////////////////////
var roll_date_a = ee.List(exceed_a.aggregate_array("system:time_start"));
var roll_date_b = ee.List(exceed_b.aggregate_array("system:time_start"));

//print (roll_date)

//either rolling mean (WARNING can exceed mem limit)
var roll2_a = ee.ImageCollection(roll_date_a
  .map(function(y) {
    var start = ee.Date(y);
    var end = start.advance(-2, 'month');
    return exceed_a.filterDate(start, end)
    .reduce(ee.Reducer.mean())
    .set('system:time_start',ee.Date(y));
  }));
  
var roll2_b = ee.ImageCollection(roll_date_b
  .map(function(y) {
    var start = ee.Date(y);
    var end = start.advance(-2, 'month');
    return exceed_b.filterDate(start, end)
    .reduce(ee.Reducer.mean())
    .set('system:time_start',ee.Date(y));
  }));

//or just the mean  
var roll1_a = exceed_a.mean();
var roll1_b = exceed_b.mean();

//choose which
var roll_a = roll1_a;
var roll_b = roll1_b;
var rolls = roll_a.addBands(roll_b);


//export results
Export.image.toAsset({
  image: rolls,
  description: 'emma_output',
  assetId: 'emma_output',
  scale: 250,
  region: aoi,
  pyramidingPolicy: {
  '.default': 'mean'
  }
});

//export results
Export.image.toDrive({
  image: roll_a,
  description: 'exceed_above',
  scale: 250,
  region: aoi,
});

Export.image.toDrive({
  image: roll_b,
  description: 'exceed_below',
  scale: 250,
  region: aoi,
});

///export ndvi data

batch.Download.ImageCollection.toAsset(ndvi_percent, 'ndvi_percent', 
      {region: aoi,
       scale: 250})

//Map.addLayer(roll_a)


