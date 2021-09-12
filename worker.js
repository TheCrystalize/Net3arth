/*helper functions*/
function compDiv(z, c) {
  const s = 1 / (c.re * c.re + c.im * c.im);
  return {
    re: (z.re * c.re + z.im * c.im) * s,
    im: (z.im * c.re - z.re * c.im) * s
  }
}

function compAdd(z, c) {
  return {
    re: z.re + c.re,
    im: z.im + c.im
  }
}

function compMult(z, c) {
  return {
    re: z.re * c.re - z.im * c.im,
    im: z.re * c.im + z.im * c.re
  }
}

function compMultScalar(z, s) {
  return {
    re: z.re * s,
    im: z.im * s
  }
}
/*Transforms*/
function mobius(a, b, c, d) {
  return function(z) {
    let ans = compDiv(compAdd(compMult(a, z), b), compAdd(compMult(c, z), d));
    return {
      ...z,
      ...ans
    }
  }
}

function scale(s) {
  return function(z) {
    let ans = compMultScalar(z, s);
    return {
      ...z,
      ...ans
    }
  }
}

function blurCircle(z){
  let a = Math.random() * Math.PI * 2,
      r = Math.sqrt(Math.random());
  let ans = {
    re: Math.cos(a) * r,
    im: Math.sin(a) * r
  }
  return {
    ...z,
    ...ans
  }
}
/*IFS stuff*/
function switchStuff(stuff, val) {
  //console.log('switch');
  //console.log(stuff);

  let total = 0;
  for(let i=0;i<stuff.length;i++){
    total+=stuff[i][0];
  }
  let rand = Math.random() * total;
  let at = 0;
  for(let i = 0; i < stuff.length;i++) {
    at+=stuff[i][0];
    if(rand < at){
      //console.log(`chose ${i}`);
      return loopStuff(stuff[i][1], val);
    }
  }
  //console.log(`didn't switch!`);
}

function loopStuff(stuff, val) {
  //console.log(stuff);

  if(typeof stuff[0] === 'string') {
    //console.log(`do ${stuff[0]}`);
    switch (stuff[0]) {
      case ("mobius"):
        return mobius(...stuff[1])(val);
      case ("scale"):
        return scale(...stuff[1])(val);
      case ("blurCircle"):
        return blurCircle(val);
    }
    throw (`${stuff[0]} not supported`);
  }

  if(!stuff[0]){return val;}

  if(typeof stuff[0][0] === 'number')
    return switchStuff(stuff, val);

  for(let i = 0; i < stuff.length; i++) {
    //console.log('loop');
    //console.log(stuff[i]);
    val = loopStuff(stuff[i], val);
  }
  return val;
}

let pointer;
let stuffToDo;
let stepsPerFrame;
let WIDTH;
let HEIGHT;
let ID;

function run() {
  let buffer = Array(WIDTH * HEIGHT).fill([0, 0, 0]);
  for(let i = 0; i < stepsPerFrame; i++) {
    //console.log(`change pointer: ${stuffToDo.main}`);
    pointer = loopStuff(stuffToDo.main, pointer);

      //console.log(`do post: ${stuffToDo.post}`);
    let val = loopStuff(stuffToDo.camera, pointer);

    if(val.re + 0.5 > 0 && val.re + 0.5 < 1 && val.im + 0.5 > 0 && val.im + 0.5 < 1) {
      let index = ((val.re + 0.5) * WIDTH >> 0) + ((val.im + 0.5) * WIDTH >> 0) * HEIGHT;
      let t = buffer[index];
      buffer[index] = [t[0] + val.red, t[1] + val.green, t[2] + val.blue];
    }
  }

  postMessage([ID,buffer]);

  stepsPerFrame = Math.min(stepsPerFrame*2, 1e7);
}

function initialize(id, job, spf, width, height){
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
}

self.onmessage = function(msg) {
  switch(msg.data[0]){
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