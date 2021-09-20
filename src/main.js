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
const canvas = document.getElementById("canvas");
canvas.style.backgroundColor = "black";
canvas.style.minWidth = WIDTH + 'px';
canvas.style.maxWidth = WIDTH + 'px';
canvas.width = WIDTH;
canvas.style.minHeight = HEIGHT + 'px';
canvas.style.maxHeight = HEIGHT + 'px';
canvas.height = HEIGHT;

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

let offscreenCanvas = canvas.transferControlToOffscreen();

let rendererThread = new Worker('src/rendererWorker.js');
rendererThread.postMessage({
  canvas: offscreenCanvas
}, [offscreenCanvas]);

function updateImage(msg) {
  rendererThread.postMessage(msg.data, [msg.data]);
}

function refreshRender() {
  consolelog("running...", "limegreen");
  for(let i = 0; i < threads.length; i++) {
    threads[i].terminate();
  }

  THREADS = parseInt(threadsUI.value);
  WIDTH = parseInt(widthUI.value);
  HEIGHT = parseInt(heightUI.value);

  canvas.style.minWidth = WIDTH + 'px';
  canvas.style.maxWidth = WIDTH + 'px';
  canvas.style.minHeight = HEIGHT + 'px';
  canvas.style.maxHeight = HEIGHT + 'px';

  rendererThread.postMessage({
    width: WIDTH,
    height: HEIGHT,
    stuffToDo: stuffToDo
  });

  let spc = STEPS_PER_CALL;

  for(let i = 0; i < THREADS; i++) {
    threads[i] = new Worker('src/worker.js');
    threads[i].postMessage(["start", [i, stuffToDo, spc, WIDTH, HEIGHT]]);
    threads[i].onmessage = updateImage;
    spc *= 1.3;
    threads[i].postMessage(["data"]);
  }
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
