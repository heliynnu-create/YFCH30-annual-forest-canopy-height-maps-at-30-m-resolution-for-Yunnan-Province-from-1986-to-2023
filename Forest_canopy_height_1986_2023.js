/**
 * @description: Generates continuous 30m forest canopy height for Yunnan (1986-2023).
 * @dependencies: GEDI L4A, Landsat Collection 2, SRTM.
 * @output: ImageCollection of canopy height.
 * @author: [Li He]
 * @date: 2026-10-19
 * @code link: https://code.earthengine.google.com/7dfeefa63e76003998bef31f9498dd9e
 */
forestage = forestage.updateMask(Attribute.lte(20));
forestage = forestage.where(forestage.eq(1986),0);
// var roi = fieldBuff;
Map.centerObject(roi,8);
var empty = ee.Image().toByte();
var outline = empty.paint({
 featureCollection:roi,  // 筛选的colletion
 color:0, //颜色透明
 width:2  //边界宽度
});
//绘制红色边界
Map.addLayer(outline, {palette: "ff0000"}, "outline");
//==============================线性插值 =============================
// 这段代码实现的是掩膜替换
function replace_mask(img, newimg, nodata) {
    nodata   = nodata || 0;
    var mask = img.mask();
    img = img.unmask(nodata);
    img = img.where(mask.not(), newimg);
    img = img.updateMask(img.neq(nodata));
    return img;
}
// 添加时间属性
/** Interpolation not considering weights */
function addTimeBand(img) {
    /** make sure mask is consistent */
    var mask = img.mask();
    var time = img.metadata('system:time_start').rename("time").mask(mask);
    return img.addBands(time);
}
// 返回提取的QA位的单波段图像,排除无用的像素（如云、阴影、水等）
var getQABits = function(image, start, end) {
    var pattern = 0;
    for (var i = start; i <= end; i++) {
       pattern += 1 << i;
    }
    // Return a single band image of the extracted QA bits
    return image.bitwiseAnd(pattern).rightShift(start);
};
// 线性插值函数
function linearInterp(imgcol, frame, nodata){
    // 时间窗口，即前后寻找的时间
    frame  = frame  || 2;
    // 缺失值，即需要插值的像素对应的值，默认是0
    nodata = nodata || 0;

    var time   = 'system:time_start';
    imgcol = imgcol.map(addTimeBand);
    
    var maxDiff = ee.Filter.maxDifference(frame * (1000*60*60*24*30), time, null, time);
    var cond    = {leftField:time, rightField:time};
    // 插值之后的影像，按照时间降序排列
    // Images after, sorted in descending order (so closest is last).
    var f1 = ee.Filter.and(maxDiff, ee.Filter.lessThanOrEquals(cond));
    var c1 = ee.Join.saveAll({matchesKey:'after', ordering:time, ascending:false})
        .apply(imgcol, imgcol, f1);
    // 插值之前的影像，按照时间升序排列
    // Images before, sorted in ascending order (so closest is last).
    var f2 = ee.Filter.and(maxDiff, ee.Filter.greaterThanOrEquals(cond));
    var c2 = ee.Join.saveAll({matchesKey:'before', ordering:time, ascending:true})
        .apply(c1, imgcol, f2);
    
    var interpolated = ee.ImageCollection(c2.map(function(img) {
        img = ee.Image(img);
        var before = ee.ImageCollection.fromImages(ee.List(img.get('before'))).mosaic();
        var after  = ee.ImageCollection.fromImages(ee.List(img.get('after'))).mosaic();
        img = img.set('before', null).set('after', null);
        // 替换空值，确认线性插值有结果
        // constrain after or before no NA values, confirm linear Interp having result
        before = replace_mask(before, after, nodata);
        after  = replace_mask(after , before, nodata);
        // 计算图像之间时间的比率
        // Compute the ratio between the image times.
        var x1 = before.select('time').double();
        var x2 = after.select('time').double();
        var now = ee.Image.constant(img.date().millis()).double();
        // 待插值影像时间减去前面的合成时间，除以后面的合成影像时间和前面的影像时间差
        var ratio = now.subtract(x1).divide(x2.subtract(x1));  // this is zero anywhere x1 = x2
        // Compute the interpolated image.
        before = before.select(0); //remove time band now;
        after  = after.select(0);
        img    = img.select(0); 
        var interp = after.subtract(before).multiply(ratio).add(before).toFloat();
        // var qc = img.mask().not().rename('qc');
        interp = replace_mask(img, interp, nodata).rename('FH');
        return interp.copyProperties(img, img.propertyNames());
    }));
    return interpolated;
}
//==============================指数计算 =============================
// Landsat8 OLI
var addVIL8 = function(img){
  var B2 = img.select("SR_B2");
  B2=B2.multiply(0.8850).add(0.0183);
  B2 = B2.rename('SR_B2');
  var B3 = img.select("SR_B3");
  B3=B3.multiply(0.9317).add(0.0123);
  B3 = B3.rename('SR_B3');
  var B4 = img.select("SR_B4");
  B4=B4.multiply(0.9372).add(0.0123);
  B4 = B4.rename('SR_B4');
  var B5 = img.select("SR_B5");
  B5=B5.multiply(0.8339).add(0.0448);
  B5 = B5.rename('SR_B5');
  var B6 = img.select("SR_B6");
  B6=B6.multiply(0.8639).add(0.0306);
  B6 = B6.rename('SR_B6');
  var B7 = img.select("SR_B7");
  B7=B7.multiply(0.9165).add(0.0116);
  B7 = B7.rename('SR_B7');
  img = ee.Image.cat([B2,B3,B4,B5,B6,B7]).float();
  var NDVI = img.normalizedDifference(['SR_B5','SR_B4']).rename('NDVI');
  var NBR = img.normalizedDifference(['SR_B5','SR_B6']).rename('NBR');
  var EVI = img.expression(
    '2.5*((NIR-Red)/(NIR+6*Red-7.5*blue+10000))',
    {
      blue: img.select('SR_B2'),    // 0.452-0.512, BLUE
      Red: img.select('SR_B4'),    // 0.636-0.673 μm, RED
      NIR: img.select('SR_B5'),    //  NIR
    }).rename('EVI');
    
  var SAVI = img.expression(
    '1.5*(NIR-Red)/(NIR+Red+0.5)',
  {
    Red: img.select('SR_B4'),    // 0.636-0.673 μm, RED
    NIR: img.select('SR_B5'),    //  NIR
  }).rename('SAVI');
  var RVI = img.select('SR_B5').divide(img.select('SR_B4')).rename("RVI");
  var NDPI = img.expression(
    '(NIR-(0.74*red + 0.26*swir1))/(NIR+(0.74*red + 0.26*swir1))',
  {
    red: img.select('SR_B4'),    // 0.636-0.673 μm, RED
    NIR: img.select('SR_B5'),    //  NIR
    swir1: img.select('SR_B6'),
  }).rename('NDPI');
  var NDGI = img.expression(
    '(0.65*NIR + 0.35*green - red)/(0.65*NIR + 0.35*green + red)',
  {
    red: img.select('SR_B4'),
    NIR: img.select('SR_B5'),
    green: img.select('SR_B3'),
  }).rename('NDGI');
  return img.addBands(NDVI)
            .addBands(NBR)
            .addBands(RVI)
            .addBands(EVI)
            .addBands(SAVI)
            .addBands(NDPI)
            .addBands(NDGI);
};
// Landsat5&7
var addVIL5 = function(img){
  var NDVI = img.normalizedDifference(['SR_B4','SR_B3']).rename('NDVI');
  var NBR = img.normalizedDifference(['SR_B4','SR_B5']).rename('NBR');
  var EVI = img.expression(
    '2.5*((NIR-Red)/(NIR+6*Red-7.5*blue+10000))',
    {
      blue: img.select('SR_B1'),    // 0.452-0.512, BLUE
      Red: img.select('SR_B3'),    // 0.636-0.673 μm, RED
      NIR: img.select('SR_B4'),    //  NIR
    }).rename('EVI');
    
  var SAVI = img.expression(
    '1.5*(NIR-Red)/(NIR+Red+0.5)',
  {
    Red: img.select('SR_B3'),    // 0.636-0.673 μm, RED
    NIR: img.select('SR_B4'),    //  NIR
  }).rename('SAVI');
  var RVI = img.select('SR_B4').divide(img.select('SR_B3')).rename("RVI");
  var NDPI = img.expression(
    '(NIR-(0.74*red + 0.26*swir1))/(NIR+(0.74*red + 0.26*swir1))',
  {
    red: img.select('SR_B3'),    // 0.636-0.673 μm, RED
    NIR: img.select('SR_B4'),    //  NIR
    swir1: img.select('SR_B5'),
  }).rename('NDPI');
  var NDGI = img.expression(
    '(0.65*NIR + 0.35*green - red)/(0.65*NIR + 0.35*green + red)',
  {
    red: img.select('SR_B3'),
    NIR: img.select('SR_B4'),
    green: img.select('SR_B2'),
  }).rename('NDGI');
  return img.addBands(NDVI)
            .addBands(NBR)
            .addBands(RVI)
            .addBands(EVI)
            .addBands(SAVI)
            .addBands(NDPI)
            .addBands(NDGI);
};
//==============================Landsat5对应的去云算法 =============================
 // 使用Landsat5 Collection 2,Level 2 QA_PIXEL波段（CFMask）去云
function cloudMaskL457(image) {
 // Bit 0 - Fill, Bit 1 - Dilated Cloud, Bit 2 - Cirrus, Bit 3 - Cloud, Bit 4 - Cloud Shadow
 var qaMask = image.select('QA_PIXEL').bitwiseAnd(parseInt('11111', 2)).eq(0);
 var saturationMask = image.select('QA_RADSAT').eq(0);
 // 将缩放系数应用于适当的波段
 var opticalBands = image.select('SR_B.').multiply(0.0000275).add(-0.2);
 var thermalBand = image.select('ST_B6').multiply(0.00341802).add(149.0);
 
 // 用缩放后的波段替换，并应用云掩膜
 return image.addBands(opticalBands, null, true)
 .addBands(thermalBand, null, true)
 .updateMask(qaMask)
 .updateMask(saturationMask);
 }
//==============================Landsat8对应的去云算法 =============================
function maskL8sr(image) {
   // Bit 0 - Fill, Bit 1 - Dilated Cloud, Bit 2 - Cirrus, Bit 3 - Cloud, Bit 4 - Cloud Shadow
   var qaMask = image.select('QA_PIXEL').bitwiseAnd(parseInt('11111', 2)).eq(0);
   var saturationMask = image.select('QA_RADSAT').eq(0);
 
   // 用缩放后的波段替换，并应用云掩膜
   var opticalBands = image.select('SR_B.').multiply(0.0000275).add(-0.2);
   var thermalBands = image.select('ST_B.*').multiply(0.00341802).add(149.0);
 
   // 用缩放后的波段替换，并应用云掩膜
   return image.addBands(opticalBands, null, true)
       .addBands(thermalBands, null, true)
       .updateMask(qaMask)
       .updateMask(saturationMask);
 }
//==============================2020年特征筛选 =============================
var dataset2020 = ee.ImageCollection('LANDSAT/LC08/C02/T1_L2')
    .filterDate('2020-01-01', '2021-10-01')
    .filterBounds(roi)
    .map(maskL8sr)
    // .map(applyScaleFactors)
    .map(addVIL8)
    .map(function(image){return image.clip(roi)})
    .median();
    // .updateMask(FC.gt(0));
print("dataset2020",dataset2020);
var nbr2020 = dataset2020.select('NBR').multiply(1000).toUint16();
var glcm2020 = nbr2020.glcmTexture({size: 1,kernel:null});
var B2_2020 = dataset2020.select("SR_B2");
var B3_2020 = dataset2020.select("SR_B3");
var B4_2020 = dataset2020.select("SR_B4");
var B5_2020 = dataset2020.select("SR_B5");
var B6_2020 = dataset2020.select("SR_B6");
var B7_2020 = dataset2020.select("SR_B7");
FH2020 = FH2020.clip(roi);
print("dataset2020",dataset2020);
//==============================标准化数值std和mean（2022年GEDI足迹点对应的统计变量） =============================
//==============================Landsat5对应的森林冠层高度数据 =============================
var years5 = ee.List.sequence(1986, 2011);
var L5_COL = ee.ImageCollection(years5.map(function(year) {
  year = ee.Number(year);
  var start = ee.Date.fromYMD(year, 1, 1);
  var end = start.advance(12, 'month');
  var tempol = ee.ImageCollection('LANDSAT/LT05/C02/T1_L2')
                // .filter(ee.Filter.calendarRange(year, year, "year")) 
                 .filterDate(start, end)
                 .map(cloudMaskL457)
                // .map(applyScaleFactors57)
                 .map(addVIL5)
                 .map(function(image){return image.clip(roi)});
  var dataset2022 = tempol.median();
  var nbr2022 = dataset2022.select('NBR').multiply(1000).toUint16();
  var glcm2022 = nbr2022.glcmTexture({size: 1,kernel:null}); 
  // 数据标准化处理
 // spectrum
  var B1 = B2_2020.subtract(dataset2022.select('SR_B1'));
  var B2 = B3_2020.subtract(dataset2022.select('SR_B2'));
  var B3 = B4_2020.subtract(dataset2022.select('SR_B3'));
  var B4 = B5_2020.subtract(dataset2022.select('SR_B4'));
  var B5 = B6_2020.subtract(dataset2022.select('SR_B5'));
  var B7 = B7_2020.subtract(dataset2022.select('SR_B7'));
  // spectrum index
  var NDVI = dataset2020.select('NDVI').subtract(dataset2022.select('NDVI'));
  var NDGI = dataset2020.select('NDGI').subtract(dataset2022.select('NDGI'));
  var NDPI = dataset2020.select('NDPI').subtract(dataset2022.select('NDPI'));
  var NBR = dataset2020.select('NBR').subtract(dataset2022.select('NBR'));
  var RVI = dataset2020.select('RVI').subtract(dataset2022.select('RVI'));
  var SAVI = dataset2020.select('SAVI').subtract(dataset2022.select('SAVI'));
  // glcm
  var asm = glcm2020.select('NBR_asm').subtract(glcm2022.select('NBR_asm'));
  var corr = glcm2020.select('NBR_corr').subtract(glcm2022.select('NBR_corr'));
  var dent = glcm2020.select('NBR_dent').subtract(glcm2022.select('NBR_dent'));
  var diss = glcm2020.select('NBR_diss').subtract(glcm2022.select('NBR_diss'));
  var dvar = glcm2020.select('NBR_dvar').subtract(glcm2022.select('NBR_dvar'));
  var ent = glcm2020.select('NBR_ent').subtract(glcm2022.select('NBR_ent'));
  var idm = glcm2020.select('NBR_idm').subtract(glcm2022.select('NBR_idm'));
  var imcorr1 = glcm2020.select('NBR_imcorr1').subtract(glcm2022.select('NBR_imcorr1'));
  var inertia = glcm2020.select('NBR_inertia').subtract(glcm2022.select('NBR_inertia'));
  var prom = glcm2020.select('NBR_prom').subtract(glcm2022.select('NBR_prom'));
  var savg = glcm2020.select('NBR_savg').subtract(glcm2022.select('NBR_savg'));
  var shade = glcm2020.select('NBR_shade').subtract(glcm2022.select('NBR_shade'));
  var svar = glcm2020.select('NBR_svar').subtract(glcm2022.select('NBR_svar'));
  var var1 = glcm2020.select('NBR_var').subtract(glcm2022.select('NBR_var'));
  // 计算差值的均值和标准差
  var img2020 = ee.Image.cat([B1, RVI, NDVI, B3, shade, B5, NBR, dvar, prom, B7, NDPI, B2, diss, B4, 
                           savg, NDGI, svar, imcorr1, asm, dent, idm, corr, ent, var1, SAVI, inertia]).float();
  var meanDict = img2020.reduceRegion({
    reducer: ee.Reducer.mean(),
    geometry: field2023, 
    tileScale: 4,
    scale: 30,
    maxPixels: 1e13
  });
  var stdDict = img2020.reduceRegion({
    reducer: ee.Reducer.stdDev(),
    geometry: field2023,
    tileScale: 4,
    scale: 30,
    maxPixels: 1e13
  });
  // 标准化处理
  var means = ee.Image.constant(meanDict.values(img2020.bandNames()));
  var stds = ee.Image.constant(stdDict.values(img2020.bandNames()));
  var imgstand = img2020.subtract(means);
  imgstand = imgstand.divide(stds);
  var FH2022 = imgstand.expression(
    "4.291+(0.26*B1)+(0.2917*RVI)-(0.2734*NDVI)-(0.1634*B3)+(0.1167*shade)-(0.0601*B5)-(0.0253*NBR)-(0.0210*dvar) + (0.0290*prom) - (0.0628*B7) - (0.0534 * NDPI) +(0.0410*B2) - (0.0086*diss) - (0.0178*B4) - (0.0291*savg) + (0.0392*NDGI) - (0.0115*svar) + (0.0097*imcorr1) - (0.0133*asm) - (0.0122*den) + (0.0097 * idm) + (0.0047*corr) + (0.0076*en) - (0.0090*var) + (0.0058*SAVI) + (0.0059*inertia)",
    {
      "B1": imgstand.select("SR_B2"),
      "B2": imgstand.select("SR_B3"),
      "B3": imgstand.select("SR_B4"),
      "B4": imgstand.select("SR_B5"),
      "B5": imgstand.select("SR_B6"),
      "B7": imgstand.select("SR_B7"),
      "NBR": imgstand.select("NBR"),
      "NDGI": imgstand.select("NDGI"),
      "NDPI": imgstand.select("NDPI"),
      "NDVI": imgstand.select("NDVI"),
      "RVI": imgstand.select("RVI"),
      "SAVI": imgstand.select("SAVI"),
      // texture
      "asm": imgstand.select("NBR_asm"),
      "corr": imgstand.select("NBR_corr"),
      "den": imgstand.select("NBR_dent"),
      "diss": imgstand.select("NBR_diss"),
      "dvar": imgstand.select("NBR_dvar"),
      "en": imgstand.select("NBR_ent"),
      "idm": imgstand.select("NBR_idm"),
      "imcorr1": imgstand.select("NBR_imcorr1"),
      "inertia": imgstand.select("NBR_inertia"),
      "prom": imgstand.select("NBR_prom"),
      "savg": imgstand.select("NBR_savg"),
      "shade": imgstand.select("NBR_shade"),
      "svar": imgstand.select("NBR_svar"),
      "var": imgstand.select("NBR_var")
    }
  ).rename("FH");
  FH2022 = (FH2020.subtract(FH2022)).rename('FH');
  // glcm
  FH2022 = FH2022.updateMask(FC.eq(10));
  FH2022 = FH2022.where(FH2022.gt(83), ee.Image(16.179));
  FH2022 = FH2022.set("year", year);
  FH2022 = FH2022.set("system:index", ee.String(year.toInt()));
  // var start = ee.Date.fromYMD(year, 1, 1);
  FH2022 = FH2022.set("system:time_start", start.millis()); 
  return FH2022;
}));

// print("L5_COL",L5_COL);
// landsat7
var years7 = ee.List.sequence(2012, 2013);
var L7_COL = ee.ImageCollection(years7.map(function(year) {
  year = ee.Number(year);
  var start = ee.Date.fromYMD(year, 1, 1);
  var end = start.advance(12, 'month');
  var tempol = ee.ImageCollection('LANDSAT/LE07/C02/T1_L2')
                // .filter(ee.Filter.calendarRange(year, year, "year")) 
                .filterDate(start, end)
                 .map(cloudMaskL457)
                // .map(applyScaleFactors57)
                 .map(function(image){return image.clip(roi)})
                 .map(addVIL5);
  var dataset2022 = tempol.median();
  var nbr2022 = dataset2022.select('NBR').multiply(1000).toUint16();
  var glcm2022 = nbr2022.glcmTexture({size: 1,kernel:null}); 
  // 数据标准化处理
 // spectrum
  var B1 = B2_2020.subtract(dataset2022.select('SR_B1'));
  var B2 = B3_2020.subtract(dataset2022.select('SR_B2'));
  var B3 = B4_2020.subtract(dataset2022.select('SR_B3'));
  var B4 = B5_2020.subtract(dataset2022.select('SR_B4'));
  var B5 = B6_2020.subtract(dataset2022.select('SR_B5'));
  var B7 = B7_2020.subtract(dataset2022.select('SR_B7'));
  // spectrum index
  var NDVI = dataset2020.select('NDVI').subtract(dataset2022.select('NDVI'));
  var NDGI = dataset2020.select('NDGI').subtract(dataset2022.select('NDGI'));
  var NDPI = dataset2020.select('NDPI').subtract(dataset2022.select('NDPI'));
  var NBR = dataset2020.select('NBR').subtract(dataset2022.select('NBR'));
  var RVI = dataset2020.select('RVI').subtract(dataset2022.select('RVI'));
  var SAVI = dataset2020.select('SAVI').subtract(dataset2022.select('SAVI'));
  // glcm
  var asm = glcm2020.select('NBR_asm').subtract(glcm2022.select('NBR_asm'));
  var corr = glcm2020.select('NBR_corr').subtract(glcm2022.select('NBR_corr'));
  var dent = glcm2020.select('NBR_dent').subtract(glcm2022.select('NBR_dent'));
  var diss = glcm2020.select('NBR_diss').subtract(glcm2022.select('NBR_diss'));
  var dvar = glcm2020.select('NBR_dvar').subtract(glcm2022.select('NBR_dvar'));
  var ent = glcm2020.select('NBR_ent').subtract(glcm2022.select('NBR_ent'));
  var idm = glcm2020.select('NBR_idm').subtract(glcm2022.select('NBR_idm'));
  var imcorr1 = glcm2020.select('NBR_imcorr1').subtract(glcm2022.select('NBR_imcorr1'));
  var inertia = glcm2020.select('NBR_inertia').subtract(glcm2022.select('NBR_inertia'));
  var prom = glcm2020.select('NBR_prom').subtract(glcm2022.select('NBR_prom'));
  var savg = glcm2020.select('NBR_savg').subtract(glcm2022.select('NBR_savg'));
  var shade = glcm2020.select('NBR_shade').subtract(glcm2022.select('NBR_shade'));
  var svar = glcm2020.select('NBR_svar').subtract(glcm2022.select('NBR_svar'));
  var var1 = glcm2020.select('NBR_var').subtract(glcm2022.select('NBR_var'));
  // 计算差值的均值和标准差
  var img2020 = ee.Image.cat([B1, RVI, NDVI, B3, shade, B5, NBR, dvar, prom, B7, NDPI, B2, diss, B4, 
                           savg, NDGI, svar, imcorr1, asm, dent, idm, corr, ent, var1, SAVI, inertia]).float();
  var meanDict = img2020.reduceRegion({
    reducer: ee.Reducer.mean(),
    geometry: field2023, 
    tileScale: 4,
    scale: 30,
    maxPixels: 1e13
  });
  var stdDict = img2020.reduceRegion({
    reducer: ee.Reducer.stdDev(),
    geometry: field2023,
    tileScale: 4,
    scale: 30,
    maxPixels: 1e13
  });
  // 标准化处理
  var means = ee.Image.constant(meanDict.values(img2020.bandNames()));
  var stds = ee.Image.constant(stdDict.values(img2020.bandNames()));
  var imgstand = img2020.subtract(means);
  imgstand = imgstand.divide(stds);
  var FH2022 = imgstand.expression(
    "4.291+(0.26*B1)+(0.2917*RVI)-(0.2734*NDVI)-(0.1634*B3)+(0.1167*shade)-(0.0601*B5)-(0.0253*NBR)-(0.0210*dvar) + (0.0290*prom) - (0.0628*B7) - (0.0534 * NDPI) +(0.0410*B2) - (0.0086*diss) - (0.0178*B4) - (0.0291*savg) + (0.0392*NDGI) - (0.0115*svar) + (0.0097*imcorr1) - (0.0133*asm) - (0.0122*den) + (0.0097 * idm) + (0.0047*corr) + (0.0076*en) - (0.0090*var) + (0.0058*SAVI) + (0.0059*inertia)",
    {
      "B1": imgstand.select("SR_B2"),
      "B2": imgstand.select("SR_B3"),
      "B3": imgstand.select("SR_B4"),
      "B4": imgstand.select("SR_B5"),
      "B5": imgstand.select("SR_B6"),
      "B7": imgstand.select("SR_B7"),
      "NBR": imgstand.select("NBR"),
      "NDGI": imgstand.select("NDGI"),
      "NDPI": imgstand.select("NDPI"),
      "NDVI": imgstand.select("NDVI"),
      "RVI": imgstand.select("RVI"),
      "SAVI": imgstand.select("SAVI"),
      // texture
      "asm": imgstand.select("NBR_asm"),
      "corr": imgstand.select("NBR_corr"),
      "den": imgstand.select("NBR_dent"),
      "diss": imgstand.select("NBR_diss"),
      "dvar": imgstand.select("NBR_dvar"),
      "en": imgstand.select("NBR_ent"),
      "idm": imgstand.select("NBR_idm"),
      "imcorr1": imgstand.select("NBR_imcorr1"),
      "inertia": imgstand.select("NBR_inertia"),
      "prom": imgstand.select("NBR_prom"),
      "savg": imgstand.select("NBR_savg"),
      "shade": imgstand.select("NBR_shade"),
      "svar": imgstand.select("NBR_svar"),
      "var": imgstand.select("NBR_var")
    }
  ).rename("FH");
  // glcm
  FH2022 = (FH2020.subtract(FH2022)).rename('FH');
  FH2022 = FH2022.updateMask(FC.eq(10));
  FH2022 = FH2022.where(FH2022.gt(83), ee.Image(16.179));
  FH2022 = FH2022.set("year", year);
  FH2022 = FH2022.set("system:index", ee.String(year.toInt()));
  // var start = ee.Date.fromYMD(year, 1, 1);
  FH2022 = FH2022.set("system:time_start", start.millis()); 
  return FH2022;
}));
// print("L7_COL",L7_COL);
//==============================Landsat8 OLI对应的森林冠层高度数据 =============================
var year8 = ee.List.sequence(2014, 2023);
var L8_COL = ee.ImageCollection(year8.map(function(year) {
  year = ee.Number(year);
  var start = ee.Date.fromYMD(year, 1, 1);
  var end = start.advance(12, 'month');
  var tempol = ee.ImageCollection('LANDSAT/LC08/C02/T1_L2')
                // .filter(ee.Filter.calendarRange(year, year, "year")) 
                 .filterDate(start, end)
                 .map(maskL8sr)
                // .map(applyScaleFactors)
                 .map(function(image){return image.clip(roi)})
                 .map(addVIL8);
  var dataset2022 = tempol.median();
  var nbr2022 = dataset2022.select('NBR').multiply(1000).toUint16();
  var glcm2022 = nbr2022.glcmTexture({size: 1,kernel:null});
  // 数据标准化处理
   // spectrum
  var B1 = dataset2020.select('SR_B2').subtract(dataset2022.select('SR_B2'));
  var B2 = dataset2020.select('SR_B3').subtract(dataset2022.select('SR_B3'));
  var B3 = dataset2020.select('SR_B4').subtract(dataset2022.select('SR_B4'));
  var B4 = dataset2020.select('SR_B5').subtract(dataset2022.select('SR_B5'));
  var B5 = dataset2020.select('SR_B6').subtract(dataset2022.select('SR_B6'));
  var B7 = dataset2020.select('SR_B7').subtract(dataset2022.select('SR_B7'));
  // spectrum index
  var NDVI = dataset2020.select('NDVI').subtract(dataset2022.select('NDVI'));
  var NDGI = dataset2020.select('NDGI').subtract(dataset2022.select('NDGI'));
  var NDPI = dataset2020.select('NDPI').subtract(dataset2022.select('NDPI'));
  var NBR = dataset2020.select('NBR').subtract(dataset2022.select('NBR'));
  var RVI = dataset2020.select('RVI').subtract(dataset2022.select('RVI'));
  var SAVI = dataset2020.select('SAVI').subtract(dataset2022.select('SAVI'));
  // glcm
  var asm = glcm2020.select('NBR_asm').subtract(glcm2022.select('NBR_asm'));
  var corr = glcm2020.select('NBR_corr').subtract(glcm2022.select('NBR_corr'));
  var dent = glcm2020.select('NBR_dent').subtract(glcm2022.select('NBR_dent'));
  var diss = glcm2020.select('NBR_diss').subtract(glcm2022.select('NBR_diss'));
  var dvar = glcm2020.select('NBR_dvar').subtract(glcm2022.select('NBR_dvar'));
  var ent = glcm2020.select('NBR_ent').subtract(glcm2022.select('NBR_ent'));
  var idm = glcm2020.select('NBR_idm').subtract(glcm2022.select('NBR_idm'));
  var imcorr1 = glcm2020.select('NBR_imcorr1').subtract(glcm2022.select('NBR_imcorr1'));
  var inertia = glcm2020.select('NBR_inertia').subtract(glcm2022.select('NBR_inertia'));
  var prom = glcm2020.select('NBR_prom').subtract(glcm2022.select('NBR_prom'));
  var savg = glcm2020.select('NBR_savg').subtract(glcm2022.select('NBR_savg'));
  var shade = glcm2020.select('NBR_shade').subtract(glcm2022.select('NBR_shade'));
  var svar = glcm2020.select('NBR_svar').subtract(glcm2022.select('NBR_svar'));
  var var1 = glcm2020.select('NBR_var').subtract(glcm2022.select('NBR_var'));
  // 计算差值的均值和标准差
  var img2020 = ee.Image.cat([B1, RVI, NDVI, B3, shade, B5, NBR, dvar, prom, B7, NDPI, B2, diss, B4, 
                           savg, NDGI, svar, imcorr1, asm, dent, idm, corr, ent, var1, SAVI, inertia]).float();
  var meanDict = img2020.reduceRegion({
    reducer: ee.Reducer.mean(),
    geometry: field2023, 
    tileScale: 4,
    scale: 30,
    maxPixels: 1e13
  });
  var stdDict = img2020.reduceRegion({
    reducer: ee.Reducer.stdDev(),
    geometry: field2023,
    tileScale: 4,
    scale: 30,
    maxPixels: 1e13
  });
  // 标准化处理
  var means = ee.Image.constant(meanDict.values(img2020.bandNames()));
  var stds = ee.Image.constant(stdDict.values(img2020.bandNames()));
  var imgstand = img2020.subtract(means);
  imgstand = imgstand.divide(stds);
  var FH2022 = imgstand.expression(
    "4.291+(0.26*B1)+(0.2917*RVI)-(0.2734*NDVI)-(0.1634*B3)+(0.1167*shade)-(0.0601*B5)-(0.0253*NBR)-(0.0210*dvar) + (0.0290*prom) - (0.0628*B7) - (0.0534 * NDPI) +(0.0410*B2) - (0.0086*diss) - (0.0178*B4) - (0.0291*savg) + (0.0392*NDGI) - (0.0115*svar) + (0.0097*imcorr1) - (0.0133*asm) - (0.0122*den) + (0.0097 * idm) + (0.0047*corr) + (0.0076*en) - (0.0090*var) + (0.0058*SAVI) + (0.0059*inertia)",
    {
      "B1": imgstand.select("SR_B2"),
      "B2": imgstand.select("SR_B3"),
      "B3": imgstand.select("SR_B4"),
      "B4": imgstand.select("SR_B5"),
      "B5": imgstand.select("SR_B6"),
      "B7": imgstand.select("SR_B7"),
      "NBR": imgstand.select("NBR"),
      "NDGI": imgstand.select("NDGI"),
      "NDPI": imgstand.select("NDPI"),
      "NDVI": imgstand.select("NDVI"),
      "RVI": imgstand.select("RVI"),
      "SAVI": imgstand.select("SAVI"),
      // texture
      "asm": imgstand.select("NBR_asm"),
      "corr": imgstand.select("NBR_corr"),
      "den": imgstand.select("NBR_dent"),
      "diss": imgstand.select("NBR_diss"),
      "dvar": imgstand.select("NBR_dvar"),
      "en": imgstand.select("NBR_ent"),
      "idm": imgstand.select("NBR_idm"),
      "imcorr1": imgstand.select("NBR_imcorr1"),
      "inertia": imgstand.select("NBR_inertia"),
      "prom": imgstand.select("NBR_prom"),
      "savg": imgstand.select("NBR_savg"),
      "shade": imgstand.select("NBR_shade"),
      "svar": imgstand.select("NBR_svar"),
      "var": imgstand.select("NBR_var")
    }
  ).rename("FH");
  // FH2022 = (FH2020.subtract(FH2022)).rename('FH');
  FH2022 = (FH2020.subtract(FH2022)).rename('FH');
  FH2022 = FH2022.updateMask(FH2022.gt(0));
  FH2022 = FH2022.where(FH2022.gt(83), ee.Image(16.179));
  FH2022 = FH2022.set("year", year);
  FH2022 = FH2022.set("system:index", ee.String(year.toInt()));
  // var start = ee.Date.fromYMD(year, 1, 01);
  FH2022 = FH2022.set("system:time_start", start.millis()); 
  return FH2022;
}));
// print("L8_COL",L8_COL);
// 1986~2023年森林冠层高度数据
var FHLandsatcol = ee.ImageCollection(L5_COL.merge(L7_COL).merge(L8_COL));
print("Landsatcol", FHLandsatcol);
// // print("1986", (FHLandsatcol.filterDate('1986', '1987')).mean());
// Map.addLayer(FHLandsatcol.filterDate('2022', '2023').mean(), {}, "2022");
// //==============================树高数据集线性插值 =============================
var nodata = -9999;
var frame = 8*4;
var Landsatinter = linearInterp(FHLandsatcol.select('FH'),frame,nodata);
print("Landsatinter", Landsatinter);
Map.addLayer(Landsatinter.filterDate('2022', '2023').mean(), {}, "2022inter");
//==============================保留扰动年份后的FCH时序数据 =============================
// 将发生森林扰动区域对应的森林冠层高度设置为0
var NBRsg = Landsatinter;
forestage = forestage.rename('yod');
var yearList1 = ee.List.sequence(1986, 2023);
var landcollmean = ee.ImageCollection(yearList1.map(function(y){
  var yod = forestage.clip(roi);
  // 获取目标年份的FCH数据
  var start = ee.String(ee.Number(y).toInt());
  var end = ee.String((ee.Number(y).add(1)).toInt());
  var tempol = NBRsg.filterDate(start, end).median();
  // 获取2023年与目标年份NBR数据差值
  var Temean = tempol.where(yod.eq(ee.Number(y)), 0);
  // var Temean = tempol.updateMask(yod.neq(ee.Number(y)));
  // Temean = Temean.rename('Temean');
  var start1 = ee.Date.fromYMD(y, 1, 1);
  Temean = Temean.set("year", y);
  Temean = Temean.set("system:index", start1.format('yyyy'));
  Temean = Temean.set("system:time_start", start1.millis());
  return Temean.rename('FCH');
}));
print("landcollmean",landcollmean);
//==============================森林冠层高度变化较大的像元用相近三年的中位数替代 =============================
//  获取扰动后时间序列的变化：利用当前年份、前一年和后一年三年的中值替换变化幅度较大（5m）的像元
var Landsatinter = landcollmean;
function chuli1986(year){
  // ee.Algorithms.If(year.eq(1986),)
  var start = ee.String(ee.Number(year).toInt());
  var end1 = ee.String((ee.Number(year).add(1)).toInt());
  var end2 = ee.String((ee.Number(year).add(2)).toInt());
  // var end3 = ee.String((ee.Number(year).add(3)).toInt());
  var img1986 = Landsatinter.filterDate(start, end1).mean();
  var img1987 = Landsatinter.filterDate(end1, end2).mean();
  // var img1988 = Landsatinter.filterDate(end2, end3).mean();
  var bianhualv1 = img1987.subtract(img1986).abs();
  var merge = ee.ImageCollection([img1986,img1987]);
  merge = merge.median();
  var yearcv = img1986.where(bianhualv1.gte(5), merge);
  // ee.Algorithms.If()
  yearcv = yearcv.set("year", year);
  var start1 = ee.Date.fromYMD(year, 1, 1);
  yearcv = yearcv.set("system:index", start1.format('yyyy'));
  yearcv = yearcv.set("system:time_start", start1.millis());
  return yearcv;
}
function chuli2023(year){
  // ee.Algorithms.If(year.eq(1986),)
  var start = ee.String(ee.Number(year).toInt());
  var end1 = ee.String((ee.Number(year).add(1)).toInt());
  var end2 = ee.String((ee.Number(year).subtract(1)).toInt());
  // var end3 = ee.String((ee.Number(y).add(3)).toInt());
  var img2023 = Landsatinter.filterDate(start, end1).mean();
  var img2022 = Landsatinter.filterDate(end2, end1).mean();
  var bianhualv1 = img2023.subtract(img2022).abs();
  var merge = ee.ImageCollection([img2022,img2023]);
  // var merge = img2023.addBands(img2022);
  merge = merge.median();
  var yearcv = img2023.where(bianhualv1.gte(5), merge);
  yearcv = yearcv.set("year", year);
  var start1 = ee.Date.fromYMD(year, 1, 1);
  yearcv = yearcv.set("system:index", start1.format('yyyy'));
  yearcv = yearcv.set("system:time_start", start1.millis());
  return yearcv;
}
function chuli(year){
  // ee.Algorithms.If(year.eq(1986),)
  var start = ee.String(ee.Number(year).toInt());
  var end1 = ee.String((ee.Number(year).add(1)).toInt());
  var end3 = ee.String((ee.Number(year).add(2)).toInt());
  var end2 = ee.String((ee.Number(year).subtract(1)).toInt());
  var img = Landsatinter.filterDate(start, end1).mean();
  var before = Landsatinter.filterDate(end2, start).mean();
  var after = Landsatinter.filterDate(end1, end3).mean();
  var bhlBefore = (img.subtract(before)).abs();
  var bhlAfter = (img.subtract(after)).abs();
  var merge = ee.ImageCollection([img,before,after]);
  // var merge = img.addBands(before).addBands(after);
  merge = merge.median();
  var yearcv = img.where((bhlBefore.gte(5)).or(bhlAfter.gte(5)), merge);
  yearcv = yearcv.set("year", year);
  var start1 = ee.Date.fromYMD(year, 1, 1);
  yearcv = yearcv.set("system:index", start1.format('yyyy'));
  yearcv = yearcv.set("system:time_start", start1.millis());
  return yearcv;
}
var inter1986 = [1986].map(chuli1986);
var inter2023 = [2023].map(chuli2023);
var yearList = ee.List.sequence(1987, 2022);
var internew = (yearList.map(chuli));
internew = internew.add(ee.Image(inter2023));
internew = internew.insert(0, ee.Image(inter1986));
internew = ee.ImageCollection(internew);
// var FH2023 = internew.filterDate('2022', '2023').mean();
// FH2023.rename('FH2023');
var table = roi;
// print("FH2023",FH2023);
print("internew", ee.ImageCollection(internew));
// Map.addLayer(FH2023, {}, "2022chuli");
// //==============================SG滤波 =============================
// var oeel=require('users/OEEL/lib:loadAll');
// var millis = ee.Number(365*5).multiply(1000*60*60*24);
// // Use the same maxDiffFilter we used earlier
// var maxDiffFilter = ee.Filter.maxDifference({
//   difference: millis,
//   leftField: 'system:time_start',
//   rightField: 'system:time_start'
// });

// // Use the default distanceFunction
// var distanceFunction = function(infromedImage, estimationImage) {
//   return ee.Image.constant(
//       ee.Number(infromedImage.get('system:time_start'))
//       .subtract(
//         ee.Number(estimationImage.get('system:time_start')))
//         );
//   };
// // 把SG滤波函数设置为一次函数,即线性函数
// var order = 1;
// var smoothed = oeel.ImageCollection.SavatskyGolayFilter(
//   internew, 
//   maxDiffFilter,
//   distanceFunction,
//   order,['FCH']);
// smoothed =  smoothed.select("d_0_FCH");
// print("smoothedSG", smoothed);
var multi_band_image = internew.toBands().rename(['FCH1986', 'FCH1987', 'FCH1988', 'FCH1989', 'FCH1990', 'FCH1991','FCH1992', 
'FCH1993', 'FCH1994', 'FCH1995', 'FCH1996', 'FCH1997', 'FCH1998', 'FCH1999', 'FCH2000','FCH2001', 'FCH2002', 'FCH2003',
'FCH2004', 'FCH2005', 'FCH2006', 'FCH2007', 'FCH2008', 'FCH2009', 'FCH2010', 'FCH2011', 'FCH2012','FCH2013', 'FCH2014', 
'FCH2015', 'FCH2016', 'FCH2017', 'FCH2018', 'FCH2019', 'FCH2020', 'FCH2021', 'FCH2022', 'FCH2023']);
print("multi_band_image", multi_band_image);
Map.addLayer(multi_band_image.select('FCH1986'),{},"FCH1986");
Export.image.toAsset({
  image: multi_band_image, 
  description: "FCH1986_2023_new", 
  assetId: "projects/ee-heliynnu/assets/FCH1986_2023_new", 
  // fileNamePrefix: 'lt-gee_disturbance_map', 
  region: roi, 
  scale: 30, 
  crs: 'EPSG:4326', 
  maxPixels: 1e13
});
