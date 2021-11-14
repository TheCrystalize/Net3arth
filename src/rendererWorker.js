'use strict';
importScripts('standardLib.js');
importScripts('rendererLib.js');

let mainBuffer;
let mainZBuffer;
let img;

let canvas;
let ctx;
let queue = 0;

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

function refreshRender(width, height) {
  WIDTH = width;
  HEIGHT = height;
  canvas.width = WIDTH;
  canvas.height = HEIGHT;

  mainBuffer = [
    new Float32Array(WIDTH * HEIGHT),
    new Float32Array(WIDTH * HEIGHT),
    new Float32Array(WIDTH * HEIGHT),
  ];
  mainZBuffer = new Float64Array(WIDTH * HEIGHT);

  img = new ImageData(WIDTH, HEIGHT);
  for(let i = 3; i < WIDTH * HEIGHT * 4; i += 4) {
    img.data[i] = 255;
  }

  ctx.fillRect(0, 0, WIDTH, HEIGHT);
}

function sleep(ms) {
  return new Promise(resolve => setTimeout(resolve, ms));
}

async function updateImage(msg) {
  queue++;
  let m = [
    new Float32Array(msg.data[1]),
    new Float32Array(msg.data[2]),
    new Float32Array(msg.data[3]),
    new Float64Array(msg.data[4]),
  ];
  let id = msg.data[0];

  if(m[0].length !== mainBuffer[0].length) {
    console.error(`missmatched buffer size: ${m[0].length} != ${mainBuffer[0].length}`);
    return;
  }

  for(let i = WIDTH * HEIGHT - 1; i >= 0; i--) {
    if(m[i * 4] !== 0 || m[i * 4 + 1] !== 0 || m[i * 4 + 2] !== 0 || m[i * 4 + 3] !== 0) {
      let result = buffer({
        red: mainBuffer[0][i],
        green: mainBuffer[1][i],
        blue: mainBuffer[2][i],
        alpha: 1,
        z: mainZBuffer[i]
      }, {
        red: m[0][i],
        green: m[1][i],
        blue: m[2][i],
        alpha: 1,
        z: m[3][i]
      });
      mainBuffer[0][i] = result.red;
      mainBuffer[1][i] = result.green;
      mainBuffer[2][i] = result.blue;
      mainZBuffer[i] = result.z;
    }
  }

  m[0] = undefined;
  m[1] = undefined;
  m[2] = undefined;
  m[3] = undefined;

  // draw onto canvas
  if(queue > 1){
    console.log(`${id} | SKIP`);
    queue--;
    return;
  }
  console.log(`${id} | WAIT`);
  await sleep(100);
  console.log(`${id} | RENDER`);
  const brightest = getBrightest();

  for(let i = WIDTH * HEIGHT - 1; i >= 0; i--) {
    let shaderResult = loopStuff(stuffToDo.shader, {
      re: i % WIDTH,
      im: (i / WIDTH) >> 0,
      z: mainZBuffer[i],
      red: mainBuffer[0][i] / brightest,
      green: mainBuffer[1][i] / brightest,
      blue: mainBuffer[2][i] / brightest,
      alpha: 1,
      zBuffer: mainZBuffer,
      mainBuffer: mainBuffer,
      width: WIDTH,
      height: HEIGHT
    });
    img.data[i * 4] = shaderResult.red * 255 >> 0;
    img.data[i * 4 + 1] = shaderResult.green * 255 >> 0;
    img.data[i * 4 + 2] = shaderResult.blue * 255 >> 0;
  }
  ctx.putImageData(img, 0, 0);
  await sleep(100);
  postMessage(0);
  queue--;
}

onmessage = function(msg) {
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
  } else if(msg.data.hasOwnProperty('stuffToDo')) {
    stuffToDo = msg.data.stuffToDo;

    customFunctions = stuffToDo.customFunctions;
    loadPreCompute(stuffToDo.preCompute);
    populateFunctions(stuffToDo.shader);
  } else if(msg.data.hasOwnProperty('reset')) {
    ctx.fillRect(0, 0, WIDTH, HEIGHT);
  } else {
    updateImage(msg);
  }
}
