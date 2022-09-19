'use strict';
importScripts('standardLib.js');
importScripts('rendererLib.js');

let samples = 0;
let scl = [];
let mainBuffer;

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

  if(samples >= stepsPerFrame){
    postMessage([
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
      ]);
    postMessage({
      steps: samples
    });

    mainBuffer = [
      new Float32Array(WIDTH * HEIGHT),
      new Float32Array(WIDTH * HEIGHT),
      new Float32Array(WIDTH * HEIGHT),
      new Float64Array(WIDTH * HEIGHT)
    ];

    stepsPerFrame = Math.min(stepsPerFrame * 4, 1e9);
    samples = 0;
  }
}

function drawBody(z) {
  drawPoint = drawCamera;
  let val = loopStuff(stuffToDo.camera, z);
  drawCamera(val);
  drawPoint = drawBody;
}

function run() {
  mainBuffer = [
    new Float32Array(WIDTH * HEIGHT),
    new Float32Array(WIDTH * HEIGHT),
    new Float32Array(WIDTH * HEIGHT),
    new Float64Array(WIDTH * HEIGHT)
  ];
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

  postMessage([
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
    ]);
  postMessage({
    steps: samples
  });

  stepsPerFrame = Math.min(stepsPerFrame * 4, 1e9);

  setTimeout(run, 1);
}

function initialize(id, job, spf, width, height) {
  ID = id;
  WIDTH = width;
  HEIGHT = height;
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
      run();
      break;
    default:
      //console.log('bad request:');
      //console.log(msg.data);
  }
}
