const WIDTH = 800;
const HEIGHT = 800;

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

/**
 * Updates the canvas content
 */

const stuffToDo = {
  main: [
    [
      3,
      [1 / 3, mobius.function({
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
      })],
      [2 / 3, mobius.function({
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
      })],
      [1, mobius.function({
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
      })]
    ]
  ],
  post: [
    mobius.function({
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
    }),
    scale.function(0.4)
  ]
};

function switchStuff(stuff, val) {
  let rand = Math.random();
  for(let i = 0; i < stuff[0];) {
    if(rand < stuff[++i][0])
      return loopStuff(stuff[i][1], val);
  }
}

function loopStuff(stuff, val) {
  if(typeof stuff === 'function')
    return stuff(val);

  if(typeof stuff[0] === 'number')
    return switchStuff(stuff, val);

  for(let i = 0; i < stuff.length; i++) {
    val = loopStuff(stuff[i], val);
  }
  return val;
}

/*

{re:-0.497297383621323782,im:-0.006511070947473171},
{re: 1.000000000000000000,im:+0.000000000000000000},
{re:-1.000000000000000000,im:-0.000000000000000000},
{re: 1.437216112833956923,im:+0.018817344280739631},

{re:-0.000000000000000000,im:-0.588229835383947423},
{re: 1.000000000000000000,im:+0.000000000000000000},
{re: 1.000000000000000000,im:+0.000000000000000000},
{re: 0.000000000000000000,im:-1.700015775886789767}

{re: 1.000000000000000000,im:+0.000000000000000000},
{re: 0.000000000000000000,im:-0.588229835383947423},
{re: 0.000000000000000000,im:-1.700015775886789767},
{re: 1.000000000000000000,im:+0.000000000000000000}

 */

let buffer = Array(WIDTH * HEIGHT).fill([0, 0, 0]);
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

let pointer = {
  re: 0.0001,
  im: 0.0001,
  red: 1,
  green: 1,
  blue: 1
};

function draw() {
  const start = Date.now();
  do {
    for(let i = 0; i < 1000; i++) {

      /* STUFF */
      pointer = loopStuff(stuffToDo.main, pointer);

      /* POST-STUFF */
      let val = loopStuff(stuffToDo.post, pointer);

      /* DRAW STUFF */
      if(val.re + 0.5 > 0 && val.re + 0.5 < 1 && val.im + 0.5 > 0 && val.im + 0.5 < 1) {
        let index = ((val.re + 0.5) * WIDTH >> 0) + ((val.im + 0.5) * WIDTH >> 0) * HEIGHT;
        let t = buffer[index];
        buffer[index] = [t[0] + val.red, t[1] + val.green, t[2] + val.blue];
      }
    }
  } while(Date.now() - start < 20);

  const brightest = getBrightest();

  for(let b, i = 0; i < WIDTH * HEIGHT; i++) {
    b = buffer[i];
    img.data[i * 4] = Math.log(b[0]) / Math.log(brightest) * 255 >> 0;
    img.data[i * 4 + 1] = Math.log(b[1]) / Math.log(brightest) * 255 >> 0;
    img.data[i * 4 + 2] = Math.log(b[2]) / Math.log(brightest) * 255 >> 0;
  }
  ctx.putImageData(img, 0, 0);

  window.requestAnimationFrame(draw);
}

draw();
