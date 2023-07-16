'use strict';
importScripts('standardLib.js');
importScripts('rendererLib.js');

let mainBuffer;
let mainZBuffer;
let img;

let canvas;
let ctx;
let start = true;

let threads = 1;
let seenThreads = 0;

render = false;
frame = 0;
frames = 1;
let code;

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

  start = true;

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

async function updateImage(msg) {
  console.log(`Shader | rendering...`);
  mainBuffer = msg.data.mainBuffer;
  mainZBuffer = msg.data.mainZBuffer;
  for(let i = WIDTH * HEIGHT - 1; i >= 0; i--) {
    let shaderResult = loopStuff(stuffToDo.shader, {
      re: i % WIDTH,
      im: (i / WIDTH) >> 0,
      z: mainZBuffer[i],
      red: mainBuffer[0][i],
      green: mainBuffer[1][i],
      blue: mainBuffer[2][i],
      alpha: 1,
      zBuffer: mainZBuffer,
      mainBuffer: mainBuffer,
      width: WIDTH,
      height: HEIGHT
    });
    img.data[i * 4] = shaderResult.red * 255 >> 0;
    img.data[i * 4 + 1] = shaderResult.green * 255 >> 0;
    img.data[i * 4 + 2] = shaderResult.blue * 255 >> 0;
    img.data[i * 4 + 3] = shaderResult.alpha * 255 >> 0;
  }
  ctx.putImageData(img, 0, 0);
  console.log(`Shader | RENDERED`);
  await sleep(100);
  postMessage(0);

  if(render) {
    frame++;
    stuffToDo = parseEverything(code);
    loadPreCompute(stuffToDo.preCompute);
    populateFunctions(stuffToDo.shader);
  }
}

let willRender = false;

onmessage = async function(msg) {
  if(msg.data.hasOwnProperty('canvas')) {
    canvas = msg.data.canvas;
    ctx = canvas.getContext("2d", {
      alpha: false
    });

    ctx.imageSmoothingQuality = "high";
  } else if(msg.data.hasOwnProperty('width')) {
    render = msg.data.render;
    frames = msg.data.frames;
    threads = msg.data.threads;
    code = msg.data.code;
    seenThreads = 0;
    refreshRender(msg.data.width, msg.data.height);
    stuffToDo = msg.data.stuffToDo;

    customFunctions = stuffToDo.customFunctions;
    loadPreCompute(stuffToDo.preCompute);
    populateFunctions(stuffToDo.shader);
    await resolvePromises();
  } else if(msg.data.hasOwnProperty('stuffToDo')) {
    render = false;
    threads = msg.data.threads;
    seenThreads = 0;
    stuffToDo = msg.data.stuffToDo;

    customFunctions = stuffToDo.customFunctions;
    loadPreCompute(stuffToDo.preCompute);
    populateFunctions(stuffToDo.shader);
    await resolvePromises();
  } else if(msg.data.hasOwnProperty('reset')) {
    ctx.fillRect(0, 0, WIDTH, HEIGHT);
  } else {
    if(!render && willRender) return;
    willRender = true;
    await sleep(100);
    updateImage(msg);
    willRender = false;
  }
}
