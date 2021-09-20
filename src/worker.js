'use strict';
importScripts('standardLib.js');
/*IFS stuff*/
function switchStuff(stuff, val) {
  //console.log('switch');
  //console.log(stuff);

  let total = 0;
  for(let i = 0; i < stuff.length; i++) {
    total += stuff[i][0];
  }
  let rand = Math.random() * total;
  let at = 0;
  for(let i = 0; i < stuff.length; i++) {
    at += stuff[i][0];
    if(rand < at) {
      //console.log(`chose ${i}`);
      return loopStuff(stuff[i][1], val);
    }
  }
  //console.log(`didn't switch!`);
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

  if(typeof stuff[0][0] === 'number')
    return switchStuff(stuff, val);

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

  for(let i = 0; i < stepsPerFrame; i++) {
    //console.log(`change pointer: ${stuffToDo.body}`);
    pointer = loopStuff(stuffToDo.body, pointer);

    //console.log(`do post: ${stuffToDo.post}`);
    let val = loopStuff(stuffToDo.camera, pointer);

    if(val.re + 0.5 > 0 && val.re + 0.5 < 1 && val.im + 0.5 > 0 && val.im + 0.5 < 1) {
      let index = ((val.re + 0.5) * WIDTH >> 0) + ((val.im + 0.5) * HEIGHT >> 0) * WIDTH;
      buffer[index * 3] += val.red;
      buffer[index * 3 + 1] += val.green;
      buffer[index * 3 + 2] += val.blue;
    }
  }

  postMessage(buffer.buffer, [buffer.buffer]);

  stepsPerFrame = Math.min(stepsPerFrame * 4, 1e7);

  setTimeout(run, 1);
}

function populateFunctionsSwitch(job) {
  for(let i = 0; i < job.length; i++) {
    populateFunctions(job[i][1]);
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
    populateFunctionsSwitch(job);
    return;
  }

  for(let i = 0; i < job.length; i++) {
    populateFunctions(job[i]);
  }
  return;
}

function loadCustomFunctions(functions) {
  for(let f in functions) {
    globalThis[functions[f].name] = new Function(...functions[f].params.map(a => a.name), functions[f].code);
    functions[f] = globalThis[functions[f].name];
  }
}

function initialize(id, job, spf, width, height) {
  ID = id;
  WIDTH = width;
  HEIGHT = height;
  pointer = {
    re: 0.001,
    im: 0.001,
    red: 1,
    green: 1,
    blue: 1
  };
  stuffToDo = job;
  stepsPerFrame = spf;
  customFunctions = stuffToDo.customFunctions;
  loadCustomFunctions(customFunctions);
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
