let WIDTH = 800;
let HEIGHT = 800;
let THREADS = (Math.max(navigator.hardwareConcurrency, 8) - 4) || 4;
const STEPS_PER_CALL = WIDTH * HEIGHT / 2;

let threadsUI = document.getElementById('threads');
let widthUI = document.getElementById('width');
let heightUI = document.getElementById('height');
threadsUI.value = THREADS;
widthUI.value = WIDTH;
heightUI.value = HEIGHT;

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

let buffer;
let img;

function getBrightest() {
  let brightest = 0;
  for(let i = 0; i < WIDTH * HEIGHT * 3; i++) {
    if(buffer[i] > brightest) {
      brightest = buffer[i];
    }
  }
  return brightest;
}

let threads = [];

let drewFrames;

function refreshRender() {
  consolelog("running...", "limegreen");
  for(let i = 0; i < threads.length; i++) {
    threads[i].terminate();
  }

  drewFrames = 0;

  THREADS = parseInt(threadsUI.value);
  WIDTH = parseInt(widthUI.value);
  HEIGHT = parseInt(heightUI.value);

  c.style.minWidth = WIDTH + 'px';
  c.style.maxWidth = WIDTH + 'px';
  c.width = WIDTH;
  c.style.minHeight = HEIGHT + 'px';
  c.style.maxHeight = HEIGHT + 'px';
  c.height = HEIGHT;

  buffer = new Uint32Array(WIDTH * HEIGHT * 3);
  for(let i = 0; i < buffer.length; i++) {
    buffer[i] = 0;
  }

  img = new ImageData(WIDTH, HEIGHT);
  for(let i = 3; i < WIDTH * HEIGHT * 4; i += 4) {
    img.data[i] = 255;
  }

  let spc = STEPS_PER_CALL;

  for(let i = 0; i < THREADS; i++) {
    threads[i] = new Worker('src/worker.js');
    threads[i].postMessage(["start", [i, stuffToDo, spc, WIDTH, HEIGHT]]);
    threads[i].onmessage = updateImage;
    spc *= 1.3;
    threads[i].postMessage(["data"]);
  }

  ctx.fillRect(0,0,WIDTH, HEIGHT);
}

function runCode() {
  compileButton.innerText = 'Compile';
  for(let i = 0; i < threads.length; i++) {
    threads[i].terminate();
  }
  run3arthLang(editor.getValue());
}

function stopCode() {
  compileButton.innerText = 'Compile';
  for(let i = 0; i < threads.length; i++) {
    threads[i].terminate();
  }
  compile3arthLang(editor.getValue());
}

function modMinute() {
  return ('' + (Date.now() % 60000) / 1000).padEnd(6, 0);
}

let onId = 1;

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
    img.data[i * 4] = buffer[i * 3] / brightest * 255 >> 0;
    img.data[i * 4 + 1] = buffer[i * 3 + 1] / brightest * 255 >> 0;
    img.data[i * 4 + 2] = buffer[i * 3 + 2] / brightest * 255 >> 0;
  }
  ctx.putImageData(img, 0, 0);
}
