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
    ans += loopToString(data[i], tab) + (i < data.length - 1 ? '\n' + tab + '-> ' : '');
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
  let customFunctions = customTransformsToString(data.customTransforms);
  let main = loopToString(data.main);
  let camera = loopToString(data.camera);
  return `${customFunctions}body: ${main};\n\ncamera: ${camera};`;
}

/*from*/
function getTypeOfWord(word){
  switch(true){
    case (word[0] === '/' && word[2] === '/'): return 'comment';
    case (word === '->'):
    case (word === '"'):
    case (word === "'"):
    case (word === '{'):
    case (word === '}'):
    case (word === '('):
    case (word === ')'):
    case (word === '['):
    case (word === ']'):
    case (word === ':'):
    case (word === ';'):
    case (word === ','):
    case (word === '+'):
    case (word === '-'):
    case (word === '='): return word;
    case (word === 'number'):
    case (word === 'boolean'):
    case (word === 'string'):
    case (word === 'complex'): return 'type';
    case (word === 'true'):
    case (word === 'false'): return 'bool';
    case (word === 'body'): return 'body';
    case (word === 'camera'): return 'camera';
    case (word === 'i'): return 'number';
    case (word[0].search(/\d|\./) >= 0): return 'number';
    default: return 'word';
  }
}

function lineToWords(line) {
  let words = [];
  let currentWord = '';
  let typeOfWord = '';
  let isCode = false;
  let bracketDepth = 0;
  let startAt = 0;
  function pushWord(word, i){
    if(word !== ''){words.push({word:word, at:startAt});}
    currentWord = '';
    startAt = i;
  }
  for(let i = 0; i < line.length; i++) {
    let char = line[i];
    switch (true) {
      case (char === '#'):
        return words;
      case (char === '{'):
      case (char === '}'):
      case (char === '('):
      case (char === ')'):
      case (char === ':'):
      case (char === ';'):
      case (char === ','):
      case (char === '+'):
      case (char === '='):
        pushWord(currentWord, i);
        pushWord(char, i+1);
        typeOfWord = '';
        break;
      case (char === '['):
      case (char === ']'):
        pushWord(currentWord, i);
        typeOfWord = '';
        break;
      case (char === '-'):
        pushWord(currentWord, i);
        currentWord = char;
        typeOfWord = 'minus';
        break;
      case (typeOfWord === 'minus' && char === '>'):
        pushWord('->', i+1);
        typeOfWord = '';
        break;
      case (char === ' '):
      case (char === '\t'):
      case (char === '\n'):
      case (char === '\r'):
        pushWord(currentWord, i+1);
        typeOfWord = '';
        break;
      case (typeOfWord === '' && char === '.' && line[i+1].search(/\d/) === 0):
        currentWord = '0.';
        typeOfWord = 'decimal';
        break;
      case ((typeOfWord === 'number' || typeOfWord === 'decimal') && char.search(/\d/) === 0):
        currentWord += char;
        break;
      case (typeOfWord === 'number' && char === '.'):
        currentWord += char;
        typeOfWord = 'decimal';
        break;
      case (typeOfWord === '' && char.search(/[a-zA-Z_]/) === 0):
        currentWord = char;
        typeOfWord = 'word';
        break;
      case (typeOfWord !== 'word' && char.search(/\d/) === 0):
        pushWord(currentWord, i);
        currentWord = char;
        typeOfWord = 'number';
        break;
      case ((typeOfWord === 'number' || typeOfWord === 'decimal') && char === 'i'):
        pushWord(currentWord+'i', i+1);
        typeOfWord = '';
        break;
      case (typeOfWord === 'word' && char.search(/\w/) === 0):
        currentWord += char;
        break;
      default:
        pushWord(currentWord, i);
        pushWord(char, i+1);
        typeOfWord = '';
    }
  }
  return words;
}

function _3arthError(word, lineNumber, line){
  return msg => {
    throw `ERROR: Expected a ${msg}. Instead got: ${word.word}\n` +
    `Line ${lineNumber}:${word.at}\n` +
    `${line}\n${new Array(word.at).fill(' ').join('')}^`;
  }
}

function parseEverything(code) {
  code = code.split('\n');
  let functions = [];
  let world = [];
  let camera = [];
  let parseState = [{is:"top"}];
  for(let i = 0; i < code.length; i++) {
    let words = lineToWords(code[i]);

    for(let j=0;j<words.length;j++){
      let word = words[j].word;
      let wordType = getTypeOfWord(word);
      let newError = _3arthError(words[j], i, code[i]);

      console.log(wordType + ': ' + word);

      switch(parseState[0].is){
        case ':':
        case "=":
          if(wordType === parseState[0].is){
            parseState.shift();
          }
          else{
            newError(parseState[0].is);
          }
        break;
        case "top":
          switch(wordType){
            case "body":
              parseState.unshift({is : "body"});
            break;
            case "camera":
              parseState.unshift({is : "camera"});
            break;
            case "word":
              parseState.unshift({is : "function", name : word});
            break;
            default:
              newError('function, "world", or "camera"');
          }
        break;
        case "function":
          if(wordType === '('){
            parseState.unshift({is: "function params", params: []});
          }
          else{
            newError('function declaration');
          }
        break;
        case "function params":
          if(wordType === 'type'){
            parseState.unshift({is: "function param", type: wordType});
            parseState.unshift({is: "function param name"});
          }
          else{
            newError('type definition');
          }
        break;
        case "function param name":
          if(wordType === 'word'){
            parseState.shift();
            parseState[0].name = word;
            parseState.unshift({is: 'function param default'});
            parseState.unshift({is: '='});
          }
          else{
            newError("parameter name");
          }
        break;
        case "function param default":
          if(wordType === 'word'){
            parseState.shift();
            parseState[0].name = word;
            parseState.unshift({is: 'function param default'});
            parseState.unshift({is: '='});
          }
          else{
            newError("parameter name");
          }
        break;
        default:
          throw [`Unhandeled state: ${parseState[0].is}`, parseState];
      }
    }
  }
}

function from3arthLang(code) {
  //try {
    parseEverything(code);
  //} catch (e) {
  //  console.error(e);
  //  alert(e);
  //}
}

getFile("klinian.3arth", r => {
  console.log(r);
  console.log('-- compile --');
  console.log(from3arthLang(r));
});
