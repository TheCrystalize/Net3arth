function sumStuff(stuff, val) {
  let total = {
    re: 0,
    im: 0,
    z: 0
  };
  for(let i = 1; i < stuff.length; i++) {
    let thisVal = loopStuff(stuff[i][1], val);
    total = {
      re: total.re + thisVal.re * stuff[i][0],
      im: total.im + thisVal.im * stuff[i][0],
      z: total.z + thisVal.z * stuff[i][0]
    };
  }
  return {
    ...val,
    ...total
  };
}

function productStuff(stuff, val) {
  let total = {
    re: 1,
    im: 1,
    z: 1,
    red: 1,
    green: 1,
    blue: 1
  };
  for(let i = 1; i < stuff.length; i++) {
    let thisVal = loopStuff(stuff[i][1], val);
    total = {
      re: total.re * thisVal.re * stuff[i][0],
      im: total.im * thisVal.im * stuff[i][0],
      z: total.z * thisVal.z * stuff[i][0],
      red: total.red * thisVal.red * stuff[i][0],
      green: total.green * thisVal.green * stuff[i][0],
      blue: total.blue * thisVal.blue * stuff[i][0]
    };
  }
  return {
    ...val,
    ...total
  };
}

function switchStuff(stuff, val) {
  let total = 0;
  for(let i = 1; i < stuff.length; i++) {
    total += stuff[i][0];
  }
  let rand = Math.random() * total;
  let at = 0;
  for(let i = 1; i < stuff.length; i++) {
    at += stuff[i][0];
    if(rand < at) {
      return loopStuff(stuff[i][1], val);
    }
  }
}

function xaosStuff(stuff, val) {
  let total = 0;
  for(let i = 1; i < stuff.length; i++) {
    total += stuff[i][1][0];
  }
  let rand = Math.random() * total;
  let at = 0;
  for(let i = 1; i < stuff.length; i++) {
    at += stuff[i][1][0];
    if(rand < at) {
      at = i;
      i = Infinity;
    }
  }

  //console.log(`start`);
  //console.log(stuff);
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
    switch (stuff[0]) {
      case 'sum':
        return sumStuff(stuff, val);
      case 'product':
        return productStuff(stuff, val);
      case 'choose':
        return switchStuff(stuff, val);
      case 'xaos':
        return xaosStuff(stuff.slice(1), val);
    }
  }

  if(!stuff[0]) {
    return val;
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

function populateFunctionsSwitch(job) {
  for(let i = 1; i < job.length; i++) {
    populateFunctions(job[i][1]);
  }
}

function populateFunctionsXaos(job) {
  for(let i = 1; i < job.length; i++) {
    populateFunctions(job[i][3]);
  }
}

const ITERATORY_FUNCTIONS = {
  'sum': populateFunctionsSwitch,
  'product': populateFunctionsSwitch,
  'choose': populateFunctionsSwitch,
  'xaos': populateFunctionsXaos,
};

function populateFunctions(job) {
  if(typeof job[0] === 'string') {
    //console.log(`JOB: ${job[0]}`);
    if(ITERATORY_FUNCTIONS.hasOwnProperty(job[0])) {
      ITERATORY_FUNCTIONS[job[0]](job);
      return;
    }
    if(BUILT_IN_TRANSFORMS.hasOwnProperty(job[0])) {
      job[0] = BUILT_IN_TRANSFORMS[job[0]](...job[1]);
      return;
    } else if(customFunctions.hasOwnProperty(job[0])) {
      job[0] = globalThis[job[0]](...job[1]);
      return;
    }
    throw (`${job[0]} not supported`);
  }

  if(!job[0]) {
    return;
  }

  for(let i = 0; i < job.length; i++) {
    populateFunctions(job[i]);
  }
  return;
}

function loadPreCompute(stuff) {
  for(let thing in stuff) {
    switch (stuff[thing].is) {
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
