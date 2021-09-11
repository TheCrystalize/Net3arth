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

function doStuff(val) {
  let rand = Math.random();



  switch (true) {
    case rand < 1 / 3:
      val = mobius.function({
        re: 0,
        im: 0
      }, {
        re: 1,
        im: 0
      }, {
        re: -1,
        im: 0
      }, {
        re: 1,
        im: 0
      })(val);
      break;
    case rand < 2 / 3:

      val = mobius.function({
        re: 0,
        im: 0
      }, {
        re: 1,
        im: 0
      }, {
        re: 1,
        im: 0
      }, {
        re: 0,
        im: 1
      })(val);
      break;
    default:

      val = mobius.function({
        re: 0,
        im: 0
      }, {
        re: 1,
        im: 0
      }, {
        re: -1,
        im: 0
      }, {
        re: 0,
        im: 0
      })(val);
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
  for(let i = 0; i < buffer.length; i++) {
    if(buffer[i][0] > brightest) {
      brightest = buffer[i][0];
    }
    if(buffer[i][1] > brightest) {
      brightest = buffer[i][1];
    }
    if(buffer[i][2] > brightest) {
      brightest = buffer[i][2];
    }
  }
  return brightest;
}

function draw() {
  let start = Date.now();
  let count = 0;
  do {
    count++;

    let val = {
      re: 0,
      im: 0,
      red: 1,
      green: 1,
      blue: 1
    };
    for(let i = 0; i < 10000; i++) {
      val = doStuff(val);
      if(val.re > -0.5 && val.re < 0.5 && val.im > -0.5 && val.im < 0.5) {
        let index = ((val.re + 0.5) * WIDTH >> 0) + ((-val.im + 0.5) * WIDTH >> 0) * HEIGHT;
        let t = buffer[index];
        buffer[index] = [t[0] + val.red, t[1] + val.green, t[2] + val.blue];
      }
    }
  } while(Date.now() - start < 12);

  let brightest = getBrightest();

  for(let i = 0; i < WIDTH * HEIGHT; i++) {
    img.data[i * 4] = buffer[i][0] * 10 // / brightest * 255;
    img.data[i * 4 + 1] = buffer[i][1] * 10 // / brightest * 255;
    img.data[i * 4 + 2] = buffer[i][2] * 10 // / brightest * 255;
  }
  ctx.putImageData(img, 0, 0);

  //ctx.fillStyle = 'red';
  //ctx.fillRect(Math.random()*WIDTH, 0, WIDTH, HEIGHT);

  window.requestAnimationFrame(draw);
}

draw();
