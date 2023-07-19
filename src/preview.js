let previewFunctions = document.getElementById('previewFunctionFunction');
let previewTemplate = document.getElementById('previewFunctionTemplate');
let previewFunctionParameters = document.getElementById('previewFunctionParameters');
let previewFunctionDescription = document.getElementById('previewFunctionDescription');

for(let i in BUILT_IN_TRANSFORMS_PARAMS) {
  let optionElement = document.createElement("option");
  optionElement.textContent = i;
  previewFunctions.children[0].appendChild(optionElement);
}

let PWIDTH = 300;
let PHEIGHT = 300;

const pcanvas = document.getElementById('previewCanvas');
pcanvas.style.backgroundColor = 'black';
pcanvas.style.minWidth = PWIDTH + 'px';
pcanvas.style.maxWidth = PWIDTH + 'px';
pcanvas.width = PWIDTH;
pcanvas.style.minHeight = PHEIGHT + 'px';
pcanvas.style.maxHeight = PHEIGHT + 'px';
pcanvas.height = PHEIGHT;

poffscreenCanvas = pcanvas.transferControlToOffscreen();
prendererThread = new Worker('src/previewWorker.js');
prendererThread.postMessage({
  canvas: poffscreenCanvas
}, [poffscreenCanvas]);

let examplePrograms = {
  normal: ['Checkered', 'Checkered Plane', '3D Teapot'],
  shader: ['Checkered Shader', 'Checkered Plane Shader', '3D Teapot Shader', '3D Teapot Orth Shader', 'Dark Gray Shader']
};

function updatePreviewData() {
  let txt = '';
  for(let i = 0; i < BUILT_IN_TRANSFORMS_PARAMS[previewFunctions.value].length; i++) {
    switch(BUILT_IN_TRANSFORMS_PARAMS[previewFunctions.value][i].type) {
      case('number'):
        txt += `<input onchange="updatePreview()" class="paramValue" type="number" step="0.1" value="${eval(BUILT_IN_TRANSFORMS_PARAMS[previewFunctions.value][i].default)}"> <span class="paramLabel">${BUILT_IN_TRANSFORMS_PARAMS[previewFunctions.value][i].name}</span><br>`;
      break;
      case('object'):
      case('array'):
        if(typeof BUILT_IN_TRANSFORMS_PARAMS[previewFunctions.value][i].default === 'string') {
          txt += `<textarea onchange="updatePreview()" class="paramValue">${BUILT_IN_TRANSFORMS_PARAMS[previewFunctions.value][i].default}</textarea> <span class="paramLabel">${BUILT_IN_TRANSFORMS_PARAMS[previewFunctions.value][i].name}</span><br>`;
        }
        else {
          txt += `<textarea onchange="updatePreview()" class="paramValue">${JSON.stringify(BUILT_IN_TRANSFORMS_PARAMS[previewFunctions.value][i].default)}</textarea> <span class="paramLabel">${BUILT_IN_TRANSFORMS_PARAMS[previewFunctions.value][i].name}</span><br>`;
        }
      break;
      case('string'):
      case('complex'):
      case('function'):
        txt += `<textarea onchange="updatePreview()" class="paramValue">${BUILT_IN_TRANSFORMS_PARAMS[previewFunctions.value][i].default}</textarea> <span class="paramLabel">${BUILT_IN_TRANSFORMS_PARAMS[previewFunctions.value][i].name}</span><br>`;
      break;
      default:
        console.log(`preview type ${BUILT_IN_TRANSFORMS_PARAMS[previewFunctions.value][i].type} unsupported`);
    }
  }
  previewFunctionDescription.innerText = BUILT_IN_TRANSFORMS_DESCRIPTIONS[previewFunctions.value];
  previewFunctionParameters.innerHTML = txt;

  let type = BUILT_IN_TRANSFORMS_PREVIEW_CODE[previewFunctions.value].type;
  if(type === 'any') {
    previewTemplate.children[0].innerHTML = '';
    for(let type in examplePrograms) {
      previewTemplate.children[0].innerHTML += examplePrograms[type].map(a=>`<option value="${a}">${a}</option>`);
    }
  }
  else {
    previewTemplate.children[0].innerHTML = examplePrograms[type].map(a=>`<option value="${a}">${a}</option>`);
  }
  previewTemplate.value = BUILT_IN_TRANSFORMS_PREVIEW_CODE[previewFunctions.value].default;

  updatePreview();
}

previewFunctions.addEventListener('change', updatePreviewData);

previewFunctionTemplate.addEventListener('change', updatePreview);

updatePreviewData();
updatePreview();

function updatePreview() {
  try {
    prendererThread.postMessage({
      width: PWIDTH,
      height: PHEIGHT,
      stuffToDo: parseEverything(getPreviewCode())
    });
  } catch (e) {
    console.error(e);
  }
}

function getPreviewParams() {
  return [...previewFunctionParameters.children].map(a => a.value).filter(a => a !== undefined);
}

function getPreviewCode() {
  let params = getPreviewParams().join(',');
  switch(previewTemplate.value) {
    case('Checkered Plane'):
return `
checker() {
return z => {
if((z.re + 1) % 1 > 0.5 ^ (z.im + 1) % 1 > 0.5) {
return {
  ...z,
  red: 1,
  green: 1,
  blue: 0.3
};
}

return {
...z,
red: 0.3,
green: 0.1,
blue: 0
}
}
}

body:
blurSquare()
-> checker()
-> tileLog(4)
-> rotate(90)
-> tileLog(4)
-> ${previewFunctions.value}(${params});

camera:
scale(1 / 5);

shader:
normalizeColors()
-> gamma(2.2);`;
    case('3D Teapot'):
return `
buffer() {
  return zBuffer();
}

body:
blurTeapot()
-> scale3D(1/6)
-> rotate3D(0, 90, 0)
-> translate3D(-0.4, 0.5, 0)
-> ${previewFunctions.value}(${params});

camera:
rotate3D(40, 0, 0)
-> rotate3D(0, 30, 0)
-> translate3D(0, 0, 3)
-> perspective3D();

shader:
productColor{
  1: sumColor{
    1: specular(-20, -60, 0.1, 0, 50);
    1: basicEnvironment(-20, -60, 1.504, 0);
  };
  1: ambientOcclusion(10, 1);
} -> gamma(2.2);`;
    case('Checkered Shader'):
return `
checker() {
return z => {
if((z.re + 1) * 4 % 1 > 0.5 ^ (z.im + 1) * 4 % 1 > 0.5) {
return {
  ...z,
  red: 1,
  green: 1,
  blue: 0.3
};
}

return {
...z,
red: 0.3,
green: 0.1,
blue: 0
}
}
}

body:
blurSquare()
-> checker();

camera:
scale(2/3);

shader:
normalizeColors()
-> gamma(2.2)
-> ${previewFunctions.value}(${params});
`;
    case('Checkered Plane Shader'):
return `
checker() {
return z => {
if((z.re + 1) % 1 > 0.5 ^ (z.im + 1) % 1 > 0.5) {
return {
...z,
red: 1,
green: 1,
blue: 0.3
};
}

return {
...z,
red: 0.3,
green: 0.1,
blue: 0
}
}
}

body:
blurSquare()
-> checker()
-> tileLog(4)
-> rotate(90)
-> tileLog(4);

camera:
scale(1 / 5);

shader:
normalizeColors()
-> gamma(2.2)
-> ${previewFunctions.value}(${params});`;
    case('3D Teapot Shader'):
      return `
buffer() {
return zBuffer();
}

body:
blurTeapot()
-> scale3D(1/6)
-> rotate3D(0, 90, 0)
-> translate3D(-0.4, 0.5, 0);

camera:
rotate3D(40, 0, 0)
-> rotate3D(0, 30, 0)
-> translate3D(0, 0, 3)
-> perspective3D();

shader:
${previewFunctions.value}(${params});`;
    case('3D Teapot Orth Shader'):
      return `
buffer() {
return zBuffer();
}

body:
blurTeapot()
-> scale3D(1/15)
-> rotate3D(0, 90, 0)
-> translate3D(0, 0.3, 0);

camera:
rotate3D(40, 0, 0)
-> rotate3D(0, 30, 0)
-> translate3D(0, 0, 3);

shader:
${previewFunctions.value}(${params});`;
    case('Dark Gray Shader'):
      return `
buffer() {
  return averageBuffer();
}

body:
blurSquare()
-> gradient([
  colorRGB(0.13, 0.13, 0.13),
  colorRGB(0.1, 0.1, 0.1),
]);

shader:
${previewFunctions.value}(${params});
`;
    default:
      return `
checker() {
  return z => {
    if((z.re + 1) * 4 % 1 > 0.5 ^ (z.im + 1) * 4 % 1 > 0.5) {
      return {
        ...z,
        red: 1,
        green: 1,
        blue: 0.3
      };
    }

    return {
      ...z,
      red: 0.3,
      green: 0.1,
      blue: 0
    }
  }
}

body:
blurSquare()
-> checker()
-> ${previewFunctions.value}(${params});

camera:
scale(2/3);

shader:
normalizeColors()
-> gamma(2.2);
`;
  }
}
