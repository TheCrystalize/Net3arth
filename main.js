const WIDTH = 800;
const HEIGHT = 800;
const STEPS_PER_FRAME = WIDTH*HEIGHT/16;

//initialize canvas
const c = document.getElementById("canvas");
const ctx = c.getContext("2d", {
  alpha: false
});
c.style.backgroundColor = "#12161b";
c.style.minWidth = WIDTH + 'px';
c.style.maxWidth = WIDTH + 'px';
c.width = WIDTH;
c.style.minHeight = HEIGHT + 'px';
c.style.maxHeight = HEIGHT + 'px';
c.height = HEIGHT;

ctx.imageSmoothingQuality = "high";

const stuffToDo = {
  main: [
    [
      3,
      [1 / 3, ["mobius", [{
        re: -0.497297383621323782,
        im: -0.006511070947473171
      }, {
        re: 1,
        im: 0
      }, {
        re: -1,
        im: 0
      }, {
        re: 1.437216112833956923,
        im: 0.018817344280739631
      }]]],
      [2 / 3, ["mobius", [{
        re: -0,
        im: -0.588229835383947423
      }, {
        re: 1,
        im: 0
      }, {
        re: 1,
        im: 0
      }, {
        re: 0,
        im: -1.700015775886789767
      }]]],
      [1, ["mobius", [{
        re: 1,
        im: 0
      }, {
        re: 0,
        im: -0.588229835383947423
      }, {
        re: 0,
        im: -1.700015775886789767
      }, {
        re: 1,
        im: 0
      }]]]
    ]
  ],
  post: [
    ["mobius", [{
      im: 1,
      re: 0
    }, {
      im: -1,
      re: 0
    }, {
      im: 1,
      re: 0
    }, {
      im: 1,
      re: 0
    }]],
    ["scale", [0.4]]
  ]
};

let buffer = Array(WIDTH * HEIGHT);
for(let i = 0; i < buffer.length; i++) {
  buffer[i] = [0, 0, 0];
}
let img = new ImageData(WIDTH, HEIGHT);
for(let i = 3; i < WIDTH * HEIGHT * 4; i += 4) {
  img.data[i] = 255;
}

function getBrightest() {
  let brightest = 0;
  for(let b, i = 0; i < WIDTH * HEIGHT; i++) {
    b = buffer[i];

    if(b[0] > brightest) {
      brightest = b[0];
    }
    if(b[1] > brightest) {
      brightest = b[1];
    }
    if(b[2] > brightest) {
      brightest = b[2];
    }
  }
  return brightest;
}

/*run*/
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

let pointer = {
  re: 0.001,
  im: 0.001,
  red: 1,
  green: 1,
  blue: 1
};

function run() {
  for(let i = 0; i < STEPS_PER_FRAME; i++) {
    pointer = loopStuff(stuffToDo.main, pointer);

    let val = loopStuff(stuffToDo.post, pointer);

    if(val.re + 0.5 > 0 && val.re + 0.5 < 1 && val.im + 0.5 > 0 && val.im + 0.5 < 1) {
      let index = ((val.re + 0.5) * WIDTH >> 0) + ((val.im + 0.5) * WIDTH >> 0) * HEIGHT;
      let t = buffer[index];
      buffer[index] = [t[0] + val.red, t[1] + val.green, t[2] + val.blue];
    }
  }
}

let lastFrame = Date.now();
let currentFrame = Date.now();

function draw() {
  currentFrame = Date.now();
  lastFrame = currentFrame;

  run();

  const brightest = getBrightest();

  for(let b, i = 0; i < WIDTH * HEIGHT; i++) {
    b = buffer[i];
    img.data[i * 4] = Math.log(b[0]) / Math.log(brightest) * 255 >> 0;
    img.data[i * 4 + 1] = Math.log(b[1]) / Math.log(brightest) * 255 >> 0;
    img.data[i * 4 + 2] = Math.log(b[2]) / Math.log(brightest) * 255 >> 0;
  }
  ctx.putImageData(img, 0, 0);

  requestAnimationFrame(draw);
}

draw();
