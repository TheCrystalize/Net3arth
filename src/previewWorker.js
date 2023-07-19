'use strict';
importScripts('standardLib.js');
importScripts('rendererLib.js');

let samples = 0;
let allSamples = 0;
let scl = [];
let mainBuffer;

let queue = 0;
let start = true;

let canvas;
let ctx;
let code;
let img;

function getBrightest() {
  let brightest = 0;
  for(let i = WIDTH * HEIGHT; i >= 0; i--) {
    if(mainBuffer[0][i] > brightest) {
      brightest = mainBuffer[0][i];
    }
    if(mainBuffer[1][i] > brightest) {
      brightest = mainBuffer[1][i];
    }
    if(mainBuffer[2][i] > brightest) {
      brightest = mainBuffer[2][i];
    }
  }
  return brightest;
}

function drawCamera(z) {
  let val = loopStuff(scl, z);
  if(val.alpha > 0 && val.re + 0.5 > 0 && val.re + 0.5 < 1 && val.im + 0.5 > 0 && val.im + 0.5 < 1) {
    samples++;
    let index = ((val.re + 0.5) * WIDTH >> 0) + ((val.im + 0.5) * HEIGHT >> 0) * WIDTH;
    let result = buffer({
      red: mainBuffer[0][index],
      green: mainBuffer[1][index],
      blue: mainBuffer[2][index],
      alpha: 1,
      z: mainBuffer[3][index]
    }, val);
    mainBuffer[0][index] = result.red;
    mainBuffer[1][index] = result.green;
    mainBuffer[2][index] = result.blue;
    mainBuffer[3][index] = result.z;
  }
}

function drawBody(z) {
  drawPoint = drawCamera;
  let val = loopStuff(stuffToDo.camera, z);
  drawCamera(val);
  drawPoint = drawBody;
}

let idle = true;
let drawNext = -1;

function run() {
  if(idle) {
    setTimeout(run, 100);
    return;
  }

  scl = [];
  if(WIDTH > HEIGHT) {
    scl = [scale2(
      HEIGHT / WIDTH,
      1
    )];
  }
  if(WIDTH < HEIGHT) {
    scl = [scale2(
      1,
      WIDTH / HEIGHT
    )];
  }

  let t = 0;
  try {
    while(samples <= stepsPerFrame && t++ < 1e5) {
      //console.log(`change pointer: ${stuffToDo.body}`);
      drawPoint = drawBody;
      pointer = loopStuff(stuffToDo.body, pointer);

      drawPoint = drawCamera;
      //console.log(`do post: ${stuffToDo.post}`);
      let val = loopStuff(stuffToDo.camera, pointer);
      val = loopStuff(scl, val);

      if(val.alpha > 0 && val.re + 0.5 > 0 && val.re + 0.5 < 1 && val.im + 0.5 > 0 && val.im + 0.5 < 1) {
        samples++;
        let index = ((val.re + 0.5) * WIDTH >> 0) + ((val.im + 0.5) * HEIGHT >> 0) * WIDTH;
        let result = buffer({
          red: mainBuffer[0][index],
          green: mainBuffer[1][index],
          blue: mainBuffer[2][index],
          alpha: 1,
          z: mainBuffer[3][index]
        }, val);
        mainBuffer[0][index] = result.red;
        mainBuffer[1][index] = result.green;
        mainBuffer[2][index] = result.blue;
        mainBuffer[3][index] = result.z;
      }
    }
  }
  catch(e) {
    idle = true;
    setTimeout(run, 100);
    return;
  }

  if(samples >= stepsPerFrame) {
    allSamples += samples;
    samples = 0;
    if(stepsPerFrame < 2e7){
      stepsPerFrame = stepsPerFrame * 4;
    }
    else{
      stepsPerFrame = 2e7;
    }
    try {
      for(let i = WIDTH * HEIGHT - 1; i >= 0; i--) {
        let shaderResult = loopStuff(stuffToDo.shader, {
          re: i % WIDTH,
          im: (i / WIDTH) >> 0,
          red: mainBuffer[0][i],
          green: mainBuffer[1][i],
          blue: mainBuffer[2][i],
          z: mainBuffer[3][i],
          alpha: 1,
          zBuffer: mainBuffer[3],
          mainBuffer: mainBuffer,
          width: WIDTH,
          height: HEIGHT
        });
        img.data[i * 4] = shaderResult.red * 255 >> 0;
        img.data[i * 4 + 1] = shaderResult.green * 255 >> 0;
        img.data[i * 4 + 2] = shaderResult.blue * 255 >> 0;
        img.data[i * 4 + 3] = shaderResult.alpha * 255 >> 0;
      }
    } catch(e) {
      idle = true;
      setTimeout(run, 100);
      return;
    }
    ctx.putImageData(img, 0, 0);
    drawNext = allSamples + stepsPerFrame - 1;
  }

  //console.log(`sample frame: ${allSamples} SL ${Math.log2(allSamples / WIDTH / HEIGHT)}`);

  if(Math.log2(allSamples / WIDTH / HEIGHT) > 6) {
    idle = true;
  }

  setTimeout(run, 4);
}

run();

function refreshRender(width, height) {
  WIDTH = width;
  HEIGHT = height;
  stepsPerFrame = WIDTH * HEIGHT;
  canvas.width = WIDTH;
  canvas.height = HEIGHT;

  start = true;
  drawNext = -1;

  allSamples = 0;
  samples = 0;

  pointer = {
    re: 0.001,
    im: 0.001,
    z: 0,
    red: 1,
    green: 1,
    blue: 1,
    alpha: 1,
  };

  mainBuffer = [
    new Float32Array(WIDTH * HEIGHT),
    new Float32Array(WIDTH * HEIGHT),
    new Float32Array(WIDTH * HEIGHT),
    new Float64Array(WIDTH * HEIGHT)
  ];

  img = new ImageData(WIDTH, HEIGHT);
  for(let i = 3; i < WIDTH * HEIGHT * 4; i += 4) {
    img.data[i] = 255;
  }

  ctx.fillRect(0, 0, WIDTH, HEIGHT);
}

onmessage = async function(msg) {
  if(msg.data.hasOwnProperty('canvas')) {
    canvas = msg.data.canvas;
    ctx = canvas.getContext("2d", {
      alpha: false
    });

    ctx.imageSmoothingQuality = "high";
  } else if(msg.data.hasOwnProperty('width')) {
    refreshRender(msg.data.width, msg.data.height);
    stuffToDo = msg.data.stuffToDo;

    customFunctions = stuffToDo.customFunctions;
    loadPreCompute(stuffToDo.preCompute);
    populateFunctions(stuffToDo.shader);
    populateFunctions(stuffToDo.body);
    populateFunctions(stuffToDo.camera);
    idle = true;
    await resolvePromises();
    idle = false;
  } else if(msg.data.hasOwnProperty('stuffToDo')) {
    refreshRender(WIDTH, HEIGHT);
    stuffToDo = msg.data.stuffToDo;

    customFunctions = stuffToDo.customFunctions;
    loadPreCompute(stuffToDo.preCompute);
    populateFunctions(stuffToDo.shader);
    populateFunctions(stuffToDo.body);
    populateFunctions(stuffToDo.camera);
    idle = true;
    await resolvePromises();
    idle = false;
  }
}
