'use strict';
importScripts('standardLib.js');
importScripts('rendererLib.js');

let mainBuffer;
let mainZBuffer;

let queue = 0;
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

  start = true;

  mainBuffer = [
    new Float32Array(WIDTH * HEIGHT),
    new Float32Array(WIDTH * HEIGHT),
    new Float32Array(WIDTH * HEIGHT),
  ];
  mainZBuffer = new Float64Array(WIDTH * HEIGHT);
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
    if(m[0][i] !== 0 || m[1][i] !== 0 || m[2][i] !== 0 || m[3][i] !== 0) {
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

  m[0] = null;
  m[1] = null;
  m[2] = null;
  m[3] = null;

  seenThreads++;
  // draw onto canvas
  if(seenThreads < threads){
    queue--;
    return;
  }
  seenThreads = 0;

  // call shader
  console.log(`Sending to shader...`);
  shaderWorker.postMessage({
    mainBuffer,
    mainZBuffer
  });

  queue--;

  if(render) {
    frame++;
    stuffToDo = parseEverything(code);
    loadPreCompute(stuffToDo.preCompute);
    populateFunctions(stuffToDo.shader);

    mainBuffer = [
      new Float32Array(WIDTH * HEIGHT),
      new Float32Array(WIDTH * HEIGHT),
      new Float32Array(WIDTH * HEIGHT),
    ];
    mainZBuffer = new Float64Array(WIDTH * HEIGHT);
  }
}

function fromShaderWorker(msg) {
  postMessage(0);
};

const shaderWorker = new Worker('shaderWorker.js');
shaderWorker.onmessage = fromShaderWorker;

onmessage = async function(msg) {
  if(msg.data.hasOwnProperty('canvas')) {
    shaderWorker.postMessage({canvas: msg.data.canvas}, [msg.data.canvas]);
  } else if(msg.data.hasOwnProperty('width')) {
    shaderWorker.postMessage(msg.data);
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
    shaderWorker.postMessage(msg.data);
    render = false;
    threads = msg.data.threads;
    seenThreads = 0;
    stuffToDo = msg.data.stuffToDo;

    customFunctions = stuffToDo.customFunctions;
    loadPreCompute(stuffToDo.preCompute);
    populateFunctions(stuffToDo.shader);
    await resolvePromises();
  } else if(msg.data.hasOwnProperty('reset')) {
    shaderWorker.postMessage(msg.data);
  } else {
    await updateImage(msg);
    await resolvePromises();
  }
}
