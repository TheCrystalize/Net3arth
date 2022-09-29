'use strict';
importScripts('standardLib.js');
importScripts('rendererLib.js');

let samples = 0;
let scl = [];
let mainBuffer;

render = false;
frame = 0;
frames = 1;
let code;

let ready = false;

let postQueue = [];

async function postData() {
  while(render && !ready) {
    await sleep(5);
  }
  ready = false;
  postMessage(...postQueue.shift());
  postMessage(...postQueue.shift());
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
/*
  if(samples >= stepsPerFrame){
    postQueue.push([[
        ID,
        mainBuffer[0].buffer,
        mainBuffer[1].buffer,
        mainBuffer[2].buffer,
        mainBuffer[3].buffer
      ],
      [
        mainBuffer[0].buffer,
        mainBuffer[1].buffer,
        mainBuffer[2].buffer,
        mainBuffer[3].buffer
      ]]);
    postQueue.push([{
      steps: samples
    }]);
    postData();

    samples = 0;

    if(render) {
      frame++;
      stuffToDo = parseEverything(code);

      loadPreCompute(stuffToDo.preCompute);
      populateFunctions(stuffToDo.body);
      populateFunctions(stuffToDo.camera);

      if(frame >= frames) {
        return;
      }
    }
    else{
      if(stepsPerFrame < 2e6){
        stepsPerFrame = stepsPerFrame * 2;
      }
      else{
        stepsPerFrame = Math.min(stepsPerFrame * 1.2, 1e9);
      }
    }

    mainBuffer = [
      new Float32Array(WIDTH * HEIGHT),
      new Float32Array(WIDTH * HEIGHT),
      new Float32Array(WIDTH * HEIGHT),
      new Float64Array(WIDTH * HEIGHT)
    ];
  }
*/
}

function drawBody(z) {
  drawPoint = drawCamera;
  let val = loopStuff(stuffToDo.camera, z);
  drawCamera(val);
  drawPoint = drawBody;
}

function run() {
  if(postQueue.length > 0) {
    setTimeout(run, 1000);
    return;
  }

  try{
    mainBuffer = [
      new Float32Array(WIDTH * HEIGHT),
      new Float32Array(WIDTH * HEIGHT),
      new Float32Array(WIDTH * HEIGHT),
      new Float64Array(WIDTH * HEIGHT)
    ];
  }
  catch(e) {
    console.log(`${ID} | no memory`);
    mainBuffer = null;
    setTimeout(run, 1000);
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

  samples = 0;

  while(samples < stepsPerFrame) {
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

  postQueue.push([[
      ID,
      mainBuffer[0].buffer,
      mainBuffer[1].buffer,
      mainBuffer[2].buffer,
      mainBuffer[3].buffer
    ],
    [
      mainBuffer[0].buffer,
      mainBuffer[1].buffer,
      mainBuffer[2].buffer,
      mainBuffer[3].buffer
    ]]);
  postQueue.push([{
    steps: samples
  }]);
  postData();

  if(render) {
    frame++;

    if(frame >= frames) {
      return;
    }

    stuffToDo = parseEverything(code);

    loadPreCompute(stuffToDo.preCompute);
    populateFunctions(stuffToDo.body);
    populateFunctions(stuffToDo.camera);
  }
  else{
    if(stepsPerFrame < 2e6){
      stepsPerFrame = stepsPerFrame * 2;
    }
    else{
      stepsPerFrame = Math.min(stepsPerFrame * 1.2, 1e9);
    }
  }

  setTimeout(run, 10);
}

function initialize(id, job, spf, width, height, _render, _frames, _code) {
  ID = id;
  WIDTH = width;
  HEIGHT = height;
  render = _render;
  frame = 0;
  frames = _frames;
  code = _code;
  pointer = {
    re: 0.001,
    im: 0.001,
    z: 0,
    red: 1,
    green: 1,
    blue: 1,
    alpha: 1,
  };
  stuffToDo = job;
  stepsPerFrame = spf;
  customFunctions = stuffToDo.customFunctions;
  loadPreCompute(stuffToDo.preCompute);
  populateFunctions(stuffToDo.body);
  populateFunctions(stuffToDo.camera);
}

self.onmessage = async function(msg) {
  switch (msg.data[0]) {
    case "start":
      initialize(...msg.data[1]);
      await resolvePromises();
      break;
    case "data":
      await resolvePromises();
      ready = true;
      run();
      break;
    case "frame":
      ready = true;
      break;
    default:
      //console.log('bad request:');
      //console.log(msg.data);
  }
}
