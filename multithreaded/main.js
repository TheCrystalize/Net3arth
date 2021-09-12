const WIDTH = 800;
const HEIGHT = 800;
const THREADS = 8;
const STEPS_PER_CALL = WIDTH;

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

let threads = [];
let activeThreads = [];

let spc = STEPS_PER_CALL;
for(let i = 0; i < THREADS; i++) {
  threads[i] = new Worker('worker.js');
  threads[i].postMessage(["start", [i, stuffToDo, STEPS_PER_CALL, WIDTH, HEIGHT]]);
  threads[i].onmessage = updateImage;
  activeThreads[i] = false;
  spc *= 2;
}

function updateImage(msg) {
  activeThreads[msg.data[0]] = false;
  let m = msg.data[1];
  if(m.length !== buffer.length) {
    return;
  }
  for(let i = 0; i < buffer.length; i++) {
    buffer[i][0] += m[i][0];
    buffer[i][1] += m[i][1];
    buffer[i][2] += m[i][2];
  }
}

let lastFrame = Date.now();
let currentFrame = Date.now();

function draw() {
  currentFrame = Date.now();
  lastFrame = currentFrame;
  const brightest = getBrightest();

  for(let b, i = 0; i < WIDTH * HEIGHT; i++) {
    b = buffer[i];
    img.data[i * 4] = Math.log(b[0]) / Math.log(brightest) * 255 >> 0;
    img.data[i * 4 + 1] = Math.log(b[1]) / Math.log(brightest) * 255 >> 0;
    img.data[i * 4 + 2] = Math.log(b[2]) / Math.log(brightest) * 255 >> 0;
  }
  ctx.putImageData(img, 0, 0);

  for(let i = 0; i < activeThreads.length; i++) {
    if(!activeThreads[i]) {
      threads[i].postMessage(["data"]);
      activeThreads[i] = true;
    }
  }

  requestAnimationFrame(draw);
}

draw();
