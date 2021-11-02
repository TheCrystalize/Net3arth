'use strict';
importScripts('standardLib.js');
importScripts('rendererLib.js');

function run() {
  let mainBuffer = new Float64Array(WIDTH * HEIGHT * 4 + 1);

  mainBuffer[WIDTH * HEIGHT * 4] = ID;

  let scl = [];
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

  let samples = 0;

  for(let i = 0; i < stepsPerFrame; i++) {
    //console.log(`change pointer: ${stuffToDo.body}`);
    pointer = loopStuff(stuffToDo.body, pointer);

    //console.log(`do post: ${stuffToDo.post}`);
    let val = loopStuff(stuffToDo.camera, pointer);
    val = loopStuff(scl, val);

    if(val.alpha > 0 && val.re + 0.5 > 0 && val.re + 0.5 < 1 && val.im + 0.5 > 0 && val.im + 0.5 < 1) {
      samples++;
      let index = ((val.re + 0.5) * WIDTH >> 0) + ((val.im + 0.5) * HEIGHT >> 0) * WIDTH;
      let result = buffer({
        red: mainBuffer[index * 4],
        green: mainBuffer[index * 4 + 1],
        blue: mainBuffer[index * 4 + 2],
        alpha: 1,
        z: mainBuffer[index * 4 + 3]
      }, val);
      mainBuffer[index * 4] = result.red;
      mainBuffer[index * 4 + 1] = result.green;
      mainBuffer[index * 4 + 2] = result.blue;
      mainBuffer[index * 4 + 3] = result.z;
    }
  }

  postMessage(mainBuffer.buffer, [mainBuffer.buffer]);
  postMessage({
    steps: samples
  });

  stepsPerFrame = Math.min(stepsPerFrame * 4, 1e7);

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

self.onmessage = function(msg) {
  switch (msg.data[0]) {
    case "start":
      initialize(...msg.data[1]);
      break;
    case "data":
      run();
      break;
    default:
      //console.log('bad request:');
      //console.log(msg.data);
  }
}
