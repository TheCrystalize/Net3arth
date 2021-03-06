let WIDTH = 1280;
let HEIGHT = 720;
let THREADS = (Math.max(navigator.hardwareConcurrency, 8) - 4) || 4;
let STEPS_PER_CALL = WIDTH * HEIGHT / 2;

const preventKeys = 'psr';

document.addEventListener('keydown', event => {
  if(event.ctrlKey && preventKeys.indexOf(event.key) >= 0) {
    event.preventDefault();
    event.stopPropagation();
  }
});

let threadsUI = document.getElementById('threads');
let widthUI = document.getElementById('width');
let heightUI = document.getElementById('height');
threadsUI.value = THREADS;
widthUI.value = WIDTH;
heightUI.value = HEIGHT;

//initialize canvas
const canvas = document.getElementById('canvas');
canvas.style.backgroundColor = 'black';
canvas.style.minWidth = WIDTH + 'px';
canvas.style.maxWidth = WIDTH + 'px';
canvas.width = WIDTH;
canvas.style.minHeight = HEIGHT + 'px';
canvas.style.maxHeight = HEIGHT + 'px';
canvas.height = HEIGHT;

let mainBuffer;
let img;

function getBrightest() {
  let brightest = 0;
  for(let i = 0; i < WIDTH * HEIGHT * 3; i++) {
    if(mainBuffer[i] > brightest) {
      brightest = mainBuffer[i];
    }
  }
  return brightest;
}

let threads = [];

let offscreenCanvas = false;
let rendererThread = false;

// offscreenCanvas isn't supported in some browsers (*cough* FireFox *cough*)
try {
  offscreenCanvas = canvas.transferControlToOffscreen();
  rendererThread = new Worker('src/rendererWorker.js');
  rendererThread.postMessage({
    canvas: offscreenCanvas
  }, [offscreenCanvas]);
  rendererThread.onmessage = incrementAnimation;
} catch (e) {
  alert(`This uses offscreenCanvas, ` +
    `which isn't supported by your browser.
    We recommend switching to Chrome, Edge, or Opera.`);
}

let loads = 0;

let samples = 0;
let loadAnimation = "/-\\|";

function incrementAnimation() {
  loads = (loads + 1) % loadAnimation.length;
  htmlConsole.children.item(htmlConsole.children.length - 1).innerText =
    `Sample Level: ${Math.floor(Math.log2(samples / WIDTH / HEIGHT) * 1000) / 1000}\nrendering ${loadAnimation[loads]}`;
}

function updateImage(msg) {
  rendererThread.postMessage([
    msg.data[0],
    msg.data[1],
    msg.data[2],
    msg.data[3],
    msg.data[4],
  ], [
    msg.data[1],
    msg.data[2],
    msg.data[3],
    msg.data[4],
  ]);
}

function workerMessage(msg) {
  if(msg.data.hasOwnProperty('steps')) {
    samples += msg.data.steps;
  } else {
    updateImage(msg);
  }
}

function refreshRender(refreshCanvas = true) {
  for(let i = 0; i < threads.length; i++) {
    threads[i].terminate();
  }

  THREADS = parseInt(threadsUI.value);

  if(refreshCanvas) {
    samples = 0;
    WIDTH = parseInt(widthUI.value);
    HEIGHT = parseInt(heightUI.value);
    STEPS_PER_CALL = WIDTH * HEIGHT / 2;

    canvas.style.minWidth = WIDTH + 'px';
    canvas.style.maxWidth = WIDTH + 'px';
    canvas.style.minHeight = HEIGHT + 'px';
    canvas.style.maxHeight = HEIGHT + 'px';

    rendererThread.postMessage({
      width: WIDTH,
      height: HEIGHT,
      stuffToDo: stuffToDo
    });
  } else {
    rendererThread.postMessage({
      stuffToDo: stuffToDo
    });
  }

  let spc = STEPS_PER_CALL;

  for(let i = 0; i < THREADS; i++) {
    threads[i] = new Worker('src/worker.js');
    threads[i].postMessage(['start', [i, stuffToDo, spc, WIDTH, HEIGHT]]);
    threads[i].onmessage = workerMessage;
    spc *= 1.3;
    threads[i].postMessage(['data']);
  }
}

let oldCode = '';

function resetCanvas() {
  loads = 0;
  rendererThread.postMessage({
    reset: true
  });
}

function runCode() {
  try {
    consoleclear();
    compileButton.innerText = 'Pause';
    for(let i = 0; i < threads.length; i++) {
      threads[i].terminate();
    }
    compile3arthLang(editor.getValue());
    let newCode = JSON.stringify([stuffToDo.customFunctions, stuffToDo.body, stuffToDo.camera]);

    if(newCode != oldCode) {
      run3arthLang(editor.getValue());
      consolelog('Running...', 'limegreen');
      oldCode = newCode;
      resetCanvas();
    } else {
      resume3arthLang(editor.getValue());
      consolelog('Resuming...', 'limegreen');
    }
    runButton.innerText = 'restart';
  } catch (e) {}
}

function stopCode() {
  try {
    for(let i = 0; i < threads.length; i++) {
      threads[i].terminate();
    }
    if(compileButton.innerText === 'Compile') {
      consoleclear();
      oldCode = '';
      resetCanvas();
      compile3arthLang(editor.getValue());
      consolelog('Finished compiling!', 'limegreen');
      runButton.innerText = 'Run';
    } else if(runButton.innerText === 'Resume') {
      consolelog('stopped.', 'limegreen');
    } else {
      consolelog('paused.', 'limegreen');
      compileButton.innerText = 'Compile';
      runButton.innerText = 'Resume';
    }
  } catch (e) {}
}

let loadDefault = true;

(new URL(window.location.href)).searchParams.forEach((x, y) =>{
  if(y === 'fractal') {
    loadDefault = false;
    getFile(`gallery/fractals/${x}/code.3arth`, r => {
      editor.setValue(r);
      editor.moveCursorToPosition({
        row: 0,
        pos: 0
      });
      editor.clearSelection();
      editor.getSession().foldAll();
    });
  }
});

if(loadDefault) {
  getFile("src/default.3arth", r => {
    let pos = editor.getCursorPosition();
    editor.setValue(r);
    editor.moveCursorToPosition(pos);
    editor.clearSelection();
  });
}
