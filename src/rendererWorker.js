'use strict';
importScripts('standardLib.js');
importScripts('rendererLib.js');

let mainBuffer;
let zBuffer;
let img;

let canvas;
let ctx;

function getBrightest() {
  let brightest = 0;
  for(let i = WIDTH * HEIGHT * 3 - 1; i >= 0; i--) {
    if(mainBuffer[i] > brightest) {
      brightest = mainBuffer[i];
    }
  }
  return brightest;
}

function refreshRender(width, height) {
  WIDTH = width;
  HEIGHT = height;
  canvas.width = WIDTH;
  canvas.height = HEIGHT;

  mainBuffer = new Uint32Array(WIDTH * HEIGHT * 3);
  zBuffer = new Float64Array(WIDTH * HEIGHT);

  img = new ImageData(WIDTH, HEIGHT);
  for(let i = 3; i < WIDTH * HEIGHT * 4; i += 4) {
    img.data[i] = 255;
  }

  ctx.fillRect(0, 0, WIDTH, HEIGHT);
}

function updateImage(msg) {
  let m = new Float64Array(msg.data);
  let id = m[WIDTH * HEIGHT * 4];

  if((m.length - 1)*0.75 !== mainBuffer.length) {
    console.error(`missmatched buffer size: ${(m.length-1)/4} != ${mainBuffer.length/3}`);
    return;
  }

  for(let i = WIDTH * HEIGHT - 1; i >= 0; i--) {
    let result = buffer({
      red: mainBuffer[i * 3],
      green: mainBuffer[i * 3 + 1],
      blue: mainBuffer[i * 3 + 2],
      alpha: 1,
      z: zBuffer[i]
    }, {
      red: m[i * 4],
      green: m[i * 4 + 1],
      blue: m[i * 4 + 2],
      alpha: 1,
      z: m[i * 4 + 3]
    });
    mainBuffer[i * 3] = result.red;
    mainBuffer[i * 3 + 1] = result.green;
    mainBuffer[i * 3 + 2] = result.blue;
    zBuffer[i] = result.z;
  }

  // draw onto canvas
  const brightest = getBrightest();

  for(let i = WIDTH * HEIGHT - 1; i >= 0; i--) {
    let shaderResult = loopStuff(stuffToDo.shader, {
      re: i % WIDTH,
      im: (i / WIDTH) >> 0,
      z: zBuffer[i],
      red: mainBuffer[i * 3] / brightest,
      green: mainBuffer[i * 3 + 1] / brightest,
      blue: mainBuffer[i * 3 + 2] / brightest,
      alpha: 1,
      zBuffer: zBuffer,
      width: WIDTH,
      height: HEIGHT
    });
    img.data[i * 4] = shaderResult.red * 255 >> 0;
    img.data[i * 4 + 1] = shaderResult.green * 255 >> 0;
    img.data[i * 4 + 2] = shaderResult.blue * 255 >> 0;
  }
  ctx.putImageData(img, 0, 0);
  postMessage(0);
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
