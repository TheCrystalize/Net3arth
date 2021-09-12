const WIDTH = 800;
const HEIGHT = 800;
const THREADS = 6;
const STEPS_PER_CALL = WIDTH * HEIGHT / 8;

let shouldDraw = 0;
let lastDraw = Date.now();

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
    [1, ["mobius", [{
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
    [1, ["mobius", [{
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
    }]]],
    [0.1,
      [
        ["blurCircle", []],
        ["scale", [0.3035586587]],
        ["mobius", [{
            re: 1,
            im: 0
          },
          {
            re: 0,
            im: 0.3035586587
          },
          {
            re: 0,
            im: -0.3035586587
          },
          {
            re: 1,
            im: 0
          }
        ]],
        ["mobius", [{
            re: 1,
            im: 0
          },
          {
            re: -1,
            im: 0
          },
          {
            re: 1,
            im: 0
          },
          {
            re: 1,
            im: 0
          }
        ]]
      ]
    ]
  ],
  camera: [
    [
      "mobius", [{
          re: 1,
          im: 0
        },
        {
          re: -1,
          im: 0
        },
        {
          re: 1,
          im: 0
        },
        {
          re: 1,
          im: 0
        }
      ],
    ],
    ["scale", [0.4]],
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

function refreshRender(){
  for(let i=0;i<threads.length;i++){
    threads[i].terminate();
  }

  let spc = STEPS_PER_CALL;

  for(let i = 0; i < THREADS; i++) {
    threads[i] = new Worker('worker.js');
    threads[i].postMessage(["start", [i, stuffToDo, STEPS_PER_CALL, WIDTH, HEIGHT]]);
    threads[i].onmessage = updateImage;
    spc *= 2;
    threads[i].postMessage(["data"]);
  }
}
refreshRender();

function updateImage(msg) {
  //console.log(`recieved ${msg.data[0]}`);
  let m = msg.data[1];
  if(m.length !== buffer.length) {
    return;
  }
  for(let i = 0; i < buffer.length; i++) {
    buffer[i][0] += m[i][0];
    buffer[i][1] += m[i][1];
    buffer[i][2] += m[i][2];
  }
  threads[msg.data[0]].postMessage(["data"]);

  let askAt = Date.now();

  shouldDraw++;
  if(shouldDraw > 5){
    draw();
  }
  else{
    setTimeout(a=>{if(lastDraw < askAt){draw();}},1000);
  }
}

function draw() {
  //console.log(`draw ${shouldDraw}`);
  shouldDraw = 0;
  lastDraw = Date.now();
  const brightest = getBrightest();

  for(let b, i = 0; i < WIDTH * HEIGHT; i++) {
    b = buffer[i];
    img.data[i * 4] = Math.log(b[0]) / Math.log(brightest) * 255 >> 0;
    img.data[i * 4 + 1] = Math.log(b[1]) / Math.log(brightest) * 255 >> 0;
    img.data[i * 4 + 2] = Math.log(b[2]) / Math.log(brightest) * 255 >> 0;
  }
  ctx.putImageData(img, 0, 0);
}
