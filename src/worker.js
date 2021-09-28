'use strict';
importScripts('standardLib.js');
/*IFS stuff*/
function switchStuff(stuff, val) {
  let total = 0;
  for(let i = 0; i < stuff.length; i++) {
    total += stuff[i][0];
  }
  let rand = Math.random() * total;
  let at = 0;
  for(let i = 0; i < stuff.length; i++) {
    at += stuff[i][0];
    if(rand < at) {
      return loopStuff(stuff[i][1], val);
    }
  }
}

function xaosStuff(stuff, val) {
  let total = 0;
  for(let i = 0; i < stuff.length; i++) {
    total += stuff[i][1][0];
  }
  let rand = Math.random() * total;
  let at = 0;
  for(let i = 0; i < stuff.length; i++) {
    at += stuff[i][1][0];
    if(rand < at) {
      at = i;
      i = Infinity;
    }
  }
  //console.log(`start`);
  while(true) {
    //console.log(`${JSON.stringify(val)} - ${at}`);
    val = loopStuff(stuff[at][3], val);
    if(stuff[at][1][1] && Math.random() < 0.5) {
      return val;
    }
    let on = at;
    at = Math.random() * stuff[on][0];
    total = 0;
    for(let i = 0; i < stuff[on][2].length; i++) {
      total += stuff[on][2][i];
      if(at <= total) {
        at = i;
        i = Infinity;
      }
    }
  }
}

function loopStuff(stuff, val) {
  //console.log(stuff);
  if(typeof stuff[0] === 'function') {
    return stuff[0](val);
  }

  if(typeof stuff[0] === 'string') {
    throw stuff;
  }

  if(!stuff[0]) {
    return val;
  }

  if(typeof stuff[0][0] === 'number') {
    switch (stuff[0].length) {
      case 2:
        return switchStuff(stuff, val);
      case 4:
        return xaosStuff(stuff, val);
    }
  }

  for(let i = 0; i < stuff.length; i++) {
    val = loopStuff(stuff[i], val);
  }
  return val;
}

let customFunctions;
let pointer;
let stuffToDo;
let stepsPerFrame;
let WIDTH;
let HEIGHT;
let ID;

function consolelog() {}

function consoleclear() {}

function run() {
  let buffer = new Float64Array(WIDTH * HEIGHT * 3 + 1);

  buffer[WIDTH * HEIGHT * 3] = ID;

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
      buffer[index * 3] += val.red * val.alpha;
      buffer[index * 3 + 1] += val.green * val.alpha;
      buffer[index * 3 + 2] += val.blue * val.alpha;
    }
  }

  postMessage(buffer.buffer, [buffer.buffer]);
  postMessage({
    steps: samples
  });

  stepsPerFrame = Math.min(stepsPerFrame * 4, 1e7);

  setTimeout(run, 1);
}

function populateFunctionsSwitch(job) {
  for(let i = 0; i < job.length; i++) {
    populateFunctions(job[i][1]);
  }
}

function populateFunctionsXaos(job) {
  for(let i = 0; i < job.length; i++) {
    populateFunctions(job[i][3]);
  }
}

function populateFunctions(job) {
  if(typeof job[0] === 'string') {
    if(BUILT_IN_TRANSFORMS.hasOwnProperty(job[0])) {
      job[0] = BUILT_IN_TRANSFORMS[job[0]](...job[1]);
      return;
    } else if(customFunctions.hasOwnProperty(job[0])) {
      job[0] = customFunctions[job[0]](...job[1]);
      return;
    }
    throw (`${job[0]} not supported`);
  }

  if(!job[0]) {
    return;
  }

  if(typeof job[0][0] === 'number') {
    switch (job[0].length) {
      case 2:
        populateFunctionsSwitch(job);
        break;
      case 4:
        populateFunctionsXaos(job);
        break;
    }
    return;
  }

  for(let i = 0; i < job.length; i++) {
    populateFunctions(job[i]);
  }
  return;
}

function loadPreCompute(stuff) {
  for(let thing in stuff) {
    switch(stuff[thing].is){
      case 'function':
        globalThis[stuff[thing].name] = new Function(...stuff[thing].params.map(a => a.name), stuff[thing].code);
        stuff[thing] = globalThis[stuff[thing].name];
        break;
      case 'const':
        globalThis[stuff[thing].name] = eval(stuff[thing].const);
        break;
    }
  }
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
