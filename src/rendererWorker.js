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


function consolelog() {}

function consoleclear() {}

let customFunctions;
let stuffToDo;
let WIDTH;
let HEIGHT;

let buffer;
let img;

let canvas;
let ctx;

function getBrightest() {
  let brightest = 0;
  for(let i = 0; i < WIDTH * HEIGHT * 3; i++) {
    if(buffer[i] > brightest) {
      brightest = buffer[i];
    }
  }
  return brightest;
}

function refreshRender(width, height) {
  WIDTH = width;
  HEIGHT = height;
  canvas.width = WIDTH;
  canvas.height = HEIGHT;

  buffer = new Uint32Array(WIDTH * HEIGHT * 3);
  for(let i = 0; i < buffer.length; i++) {
    buffer[i] = 0;
  }

  img = new ImageData(WIDTH, HEIGHT);
  for(let i = 3; i < WIDTH * HEIGHT * 4; i += 4) {
    img.data[i] = 255;
  }

  ctx.fillRect(0, 0, WIDTH, HEIGHT);
}

function updateImage(msg) {
  let m = new Float64Array(msg.data);
  let id = m[WIDTH * HEIGHT * 3];

  if(m.length - 1 !== buffer.length) {
    return;
  }

  for(let i = 0; i < buffer.length - 1; i++) {
    buffer[i] += m[i];
    buffer[i + 1] += m[i + 1];
    buffer[i + 2] += m[i + 2];
  }

  // draw onto canvas
  const brightest = getBrightest();

  for(let b, i = 0; i < WIDTH * HEIGHT; i++) {
    b = buffer[i];
    let shaderResult = loopStuff(stuffToDo.shader, {
      re: i % WIDTH,
      im: (i / WIDTH) >> 0,
      red: buffer[i * 3] / brightest,
      green: buffer[i * 3 + 1] / brightest,
      blue: buffer[i * 3 + 2] / brightest,
    });
    img.data[i * 4] = shaderResult.red * 255 >> 0;
    img.data[i * 4 + 1] = shaderResult.green * 255 >> 0;
    img.data[i * 4 + 2] = shaderResult.blue * 255 >> 0;
  }
  ctx.putImageData(img, 0, 0);
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
    loadCustomFunctions(customFunctions);
    populateFunctions(stuffToDo.shader);
  } else if(msg.data.hasOwnProperty('stuffToDo')) {
    stuffToDo = msg.data.stuffToDo;

    customFunctions = stuffToDo.customFunctions;
    loadCustomFunctions(customFunctions);
    populateFunctions(stuffToDo.shader);
  } else {
    updateImage(msg);
  }
}
