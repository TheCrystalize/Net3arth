/*to*/
function complexToString(z) {
  if(z.im === 0) {
    return z.re;
  }
  if(z.re === 0) {
    return z.im + 'i';
  }
  if(z.im < 0) {
    return `${z.re}-${-z.im}i`;
  }
  return `${z.re}+${z.im}i`;
}

function paramsToString(params) {
  let ans = '';
  for(let i = 0; i < params.length; i++) {
    switch (typeof params[i]) {
      case 'number':
        ans += params[i];
        break;
      case 'object':
        ans += complexToString(params[i]);
        break;
      default:
        throw 'bad parameter';
    }
    if(i < params.length - 1) {
      ans += ', ';
    }
  }
  return ans;
}

function switchToString(data, tab) {
  let ans = '\n' + tab + 'choose{\n';
  for(let i = 0; i < data.length; i++) {
    ans += tab + '  ' + data[i][0] + ': ' + loopToString(data[i][1], tab + '    ') + ';\n';
  }
  return ans + tab + '}';
}

function loopToString(data, tab = '') {
  if(typeof data[0] === 'string') {
    return `${data[0]}(${paramsToString(data[1])})`;
  }
  if(!data[0]) {
    return '';
  }
  if(typeof data[0][0] === 'number') {
    return switchToString(data, tab);
  }

  let ans = '\n' + tab;
  for(let i = 0; i < data.length; i++) {
    ans += loopToString(data[i], tab) + (i < data.length - 1 ? ',\n' + tab : '');
  }
  return ans;
}

function customTransformsToString(transforms) {
  let ans = '';
  for(let f in transforms) {
    ans += transforms[f].name.replace(/ /g, '_') + '(' +
      transforms[f].parameters.map(p => p.type + ' ' + p.name + ' = ' + p.default).
      join(',') + '){\n  ' +
      transforms[f].function.replace(/\n/g, '\n  ') + '\n}\n\n';
  }
  return ans;
}

function to3arthLang(data) {
    let customFunctions=customTransformsToString(data.customTransforms);
    let main=loopToString(data.main);
    let camera=loopToString(data.camera);
    return `${customFunctions}body: ${main};\n\ncamera: ${camera};`;
}

/*from*/
function from3arthLang(code) {

}
