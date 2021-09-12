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
/*IFS stuff*/
function switchStuff(stuff, val) {
  let rand = Math.random();
  for(let i = 0; i < stuff[0];) {
    if(rand < stuff[++i][0])
      return loopStuff(stuff[i][1], val);
  }
}

function loopStuff(stuff, val) {
  if(typeof stuff[0] === 'string') {
    switch (stuff[0]) {
      case ("mobius"):
        return mobius(...stuff[1])(val);
      case ("scale"):
        return scale(...stuff[1])(val);
    }
  }

  if(typeof stuff[0] === 'number')
    return switchStuff(stuff, val);

  for(let i = 0; i < stuff.length; i++) {
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
    pointer = loopStuff(stuffToDo.main, pointer);

    let val = loopStuff(stuffToDo.post, pointer);

    if(val.re + 0.5 > 0 && val.re + 0.5 < 1 && val.im + 0.5 > 0 && val.im + 0.5 < 1) {
      let index = ((val.re + 0.5) * WIDTH >> 0) + ((val.im + 0.5) * WIDTH >> 0) * HEIGHT;
      let t = buffer[index];
      buffer[index] = [t[0] + val.red, t[1] + val.green, t[2] + val.blue];
    }
  }
  postMessage([ID,buffer]);

  stepsPerFrame *= 2;
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
      console.log('bad request:');
      console.log(msg.data);
  }
}
