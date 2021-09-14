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
function getTypeOfWord(word) {
  switch (true) {
    case (word[0] === '/' && word[2] === '/'):
      return 'comment';
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
    case (word === '='):
    case (word === 'body'):
    case (word === 'camera'):
    case (word === 'choose'):
      return word;
    case (word === 'number'):
    case (word === 'boolean'):
    case (word === 'string'):
    case (word === 'complex'):
      return 'type';
    case (word === 'true'):
    case (word === 'false'):
      return 'bool';
    case (word === 'i'):
      return 'imaginary unit';
    case (word[word.length-1] === 'i' && word[0].search(/\d/) >= 0):
      return 'imaginary';
    case (word[0].search(/\d/) >= 0):
      return 'number';
    default:
      return 'word';
  }
}

function lineToWords(line) {
  let words = [];
  let currentWord = '';
  let typeOfWord = '';
  let isCode = false;
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


const errorModal = document.getElementById("errorModal");
const errorModalContent = document.getElementById("errorModalContent");

function General3arthError(word, lineNumber, line) {
  return msg => {
    errorModalContent.innerText = `ERROR: ${msg}\n` +
      `Line ${lineNumber}:${word.at}\n` +
      `${line}\n${new Array(word.at).fill(' ').join('')}^`.replace(/ /g,'\u00a0');

    errorModal.style.display = "block";

    throw `ERROR: ${msg}\n` +
      `Line ${lineNumber}:${word.at}\n` +
      `${line}\n${new Array(word.at).fill(' ').join('')}^`;
  }
}

function _3arthError(word, lineNumber, line) {
  return msg => {
    General3arthError(word, lineNumber, line)(`Expected a ${msg}. Instead got: ${word.word}`, lineNumber, line);
  }
}

function parseEverything(code) {
  code = code.split('\n');
  let parseState = [{is:"top"}];
  let customFunctions = {};
  for(let i = 0; i < code.length; i++) {
    let words = lineToWords(code[i]);

    wordLoop:
    for(let j=0;j<words.length;j++){
      let word = words[j].word;
      let wordType = getTypeOfWord(word);
      let newError = _3arthError(words[j], i, code[i]);
      let newGeneralError = General3arthError(words[j], i, code[i]);

      //console.log(customFunctions);
      //console.log(JSON.stringify(parseState));
      //console.log(parseState[0]);

      //console.log(wordType + ': ' + word);

      if(word === '/' && words[j+1].word === '/'){
        j = words.length;
        continue;
      }

      if(word === '/' && words[j+1].word === '*'){
        let err = 'ERROR: unexpected end of input; unclosed "/*"\n' +
          `Line ${i}:${words[j].at}\n` +
          `${code[i]}\n${new Array(words[j+1].at).fill(' ').join('')}^`;
        while(true){
          j++;
          while(j >= words.length){
            j = 0;
            i++;
            if(i >= code.length){
              throw err;
            }
            words = lineToWords(code[i]);
          }
          if(words[j].word === '*' && words[j+1] && words[j+1].word === '/'){
            j++;
            continue wordLoop;
          }
        }
      }

      switch(parseState[0].is){
        case ':':
        case '=':
        case '(':
        case ')':
        case '{':
        case '}':
          if(wordType === parseState[0].is){
            parseState.shift();
          }
          else{
            newError(parseState[0].is);
          }
        break;
        case 'top':
          switch(wordType){
            case 'body':
              parseState.unshift({is: 'body'});
              parseState.unshift({is: 'transform', transforms: []});
              parseState.unshift({is: ':'});
            break;
            case 'camera':
              parseState.unshift({is: 'camera'});
              parseState.unshift({is: 'transform', transforms: []});
              parseState.unshift({is: ':'});
            break;
            case 'word':
              parseState.unshift({is : 'function', name : word});
            break;
            default:
              newError('function, "world", or "camera"');
          }
        break;
        case 'function':
          if(wordType === '('){
            parseState.unshift({is: 'function params', params: []});
          }
          else{
            newError('function declaration');
          }
        break;
        case 'function params':
          if(wordType === 'type'){
            parseState.unshift({is: 'function param', type: word});
            parseState.unshift({is: word});
            parseState.unshift({is: 'function param name'});
          }
          else{
            newError('type definition');
          }
        break;
        case 'function param name':
          if(wordType === 'word'){
            parseState.shift();
            parseState[1].name = word;
            parseState.unshift({is: '='});
          }
          else{
            newError("parameter name");
          }
        break;
        case 'function param default':
          if(wordType === 'word'){
            parseState.shift();
            parseState[0].name = word;
            parseState.unshift({is: 'function param default'});
            parseState.unshift({is: '='});
          }
          else{
            newError("parameter default");
          }
        break;
        case 'number':
          switch(wordType){
            case 'number':
              if(parseState[0].hasOwnProperty('value')){
                throw newError('"," or ")"');
              }
              else{
                parseState[0].value = word;
              }
            break;
            case ',':
              var param = {
                name: parseState[1].name,
                type: parseState[1].type,
                default: parseState[0].value,
              };
              parseState.shift();
              parseState.shift();
              parseState[0].params.push(param);
            break;
            case ')':
              var param = {
                name: parseState[1].name,
                type: parseState[1].type,
                default: parseState[0].value,
              };
              parseState.shift();
              parseState.shift();
              parseState[0].params.push(param);
              parseState[1].params = parseState[0].params;
              parseState.shift();
              parseState.unshift({is: 'js'});
            break;
            default:
              newError('number');
          }
        break;
        case 'number param':
          switch(wordType){
            case '-':
              if(parseState[0].hasOwnProperty('sign')) {
                newError('number');
              }
              parseState[0].sign = '-';
              break;
            case '+':
              if(parseState[0].hasOwnProperty('sign')) {
                newError('number');
              }
              parseState[0].sign = '+';
              break;
            case 'number':
              if(parseState[0].hasOwnProperty('value')) {
                newError('"," or ")"');
              }
              if(parseState[0].hasOwnProperty('sign')){
                parseState[0].value = parseState[0].sign + word;
                delete parseState[0].sign;
                break;
              }
              parseState[0].value = word;
            break;
            case ',':
              if(!parseState[0].hasOwnProperty('value')){
                newError('number');
              }
              var val = parseFloat(parseState[0].value);
              parseState.shift();
              parseState[0].params.push(val);
              parseState[0].on++;
              if(parseState[0].on >= parseState[0].paramTypes.length){
                newError(')');
              }
              parseState.unshift({is: parseState[0].paramTypes[parseState[0].on].type + ' param'});
            break;
            case ')':
              if(!parseState[0].hasOwnProperty('value')){
                newError('number');
              }
              var val = parseFloat(parseState[0].value);
              parseState.shift();
              parseState[0].params.push(val);
              parseState[0].on++;
              if(parseState[0].on < parseState[0].paramTypes.length){
                newError(',');
              }
              parseState[1].transforms.push([parseState[0].transform, parseState[0].params]);
              parseState.shift();
              parseState.unshift({is: 'after transform'});
            break;
            default:
            newError('number');
          }
        break;
        case 'complex param':
          switch(wordType){
            case '-':
              if(parseState[0].hasOwnProperty('sign')) {
                newError('complex number');
              }
              parseState[0].sign = '-';
              break;
            case '+':
              if(parseState[0].hasOwnProperty('sign')) {
                newError('complex number');
              }
              parseState[0].sign = '+';
              break;
            case 'number':
              if(parseState[0].hasOwnProperty('imaginary')) {
                newError('",", or ")"');
              }
              if(parseState[0].hasOwnProperty('real')) {
                newError('complex pair, ",", or ")"');
              }
              if(parseState[0].hasOwnProperty('sign')){
                parseState[0].real = parseState[0].sign + word;
                delete parseState[0].sign;
                break;
              }
              parseState[0].real = word;
              break;
            case 'imaginary':
              if(parseState[0].hasOwnProperty('imaginary')) {
                newError('",", or ")"');
              }
              if(parseState[0].hasOwnProperty('sign')){
                parseState[0].imaginary = parseState[0].sign + word;
                delete parseState[0].sign;
                break;
              }
              parseState[0].imaginary = word;
            break;
            case 'imaginary unit':
              if(parseState[0].hasOwnProperty('imaginary')) {
                newError('",", or ")"');
              }
              if(parseState[0].hasOwnProperty('sign')){
                parseState[0].imaginary = parseState[0].sign + 1;
                delete parseState[0].sign;
                break;
              }
              parseState[0].imaginary = 1;
            break;
            case ',':
              var val = C(parseState[0].hasOwnProperty('real')?parseFloat(parseState[0].real):0, parseState[0].hasOwnProperty('imaginary')?parseFloat(parseState[0].imaginary):0);
              parseState.shift();
              parseState[0].params.push(val);
              parseState[0].on++;
              if(parseState[0].on >= parseState[0].paramTypes.length){
                newError(')');
              }
              parseState.unshift({is: parseState[0].paramTypes[parseState[0].on].type + ' param'});
            break;
            case ')':
              var val = C(parseState[0].hasOwnProperty('real')?parseFloat(parseState[0].real):0, parseState[0].hasOwnProperty('imaginary')?parseFloat(parseState[0].imaginary):0);
              parseState.shift();
              parseState[0].params.push(val);
              parseState[0].on++;
              if(parseState[0].on < parseState[0].paramTypes.length){
                newError(',');
              }
              parseState[1].transforms.push([parseState[0].transform, parseState[0].params]);
              parseState.shift();
              parseState.unshift({is: 'after transform'});
            break;
            default:
              newError('complex number');
          }
        break;
        case 'js':
          if(word === '{'){
            j++;
            if(j >= words.length){
              j=0;
              i++;
              if(i >= code.length){
                throw err;
              }
              words = lineToWords(code[i]);
            }
            let err = 'ERROR: unexpected end of input; unclosed "{"\n' +
              `Line ${i}:${words[j].at}\n` +
              `${code[i]}\n${new Array(words[j].at).fill(' ').join('')}^`;
            let start = [i,words[j].at];
            let bracketDepth = 1;
            while(bracketDepth > 0){
              j++;
              while(j >= words.length){
                j=0;
                i++;
                if(i >= code.length){
                  throw err;
                }
                words = lineToWords(code[i]);
              }
              if(words[j].word === '{'){
                bracketDepth++;
              }
              if(words[j].word === '}'){
                bracketDepth--;
              }
            }
            parseState.shift();

            let jsCode = '';

            if(start[0] === i){
              jsCode = code[i].slice(start[1],words[j].at);
            }
            else{
              jsCode = code[start[0]].slice(start[1]) + '\n';
              for(let l = start[0]+1; l < i; l++){
                jsCode += code[l] + '\n';
              }
              jsCode += code[i].slice(0,words[j].at);
            }

            customFunctions[parseState[0].name] = {
              name: parseState[0].name,
              params: parseState[0].params,
              code: jsCode,
            };
            parseState.shift();
          }
          else{
            newError('JavaScript function');
          }
        break;
        case 'transform':
          switch(wordType){
            case 'choose':
              parseState.unshift({is: 'choose items', items: []});
              parseState.unshift({is: '{'});
            break;
            case 'word':
              for(let transform in BUILT_IN_TRANSFORMS) {
                if(word === transform){
                  //console.log(`found ${transform}`);
                  if(BUILT_IN_TRANSFORMS[word].length === 0){
                    parseState[0].transforms.push([word,[]]);
                    parseState.unshift({is: 'after transform'});
                    parseState.unshift({is: ')'});
                    parseState.unshift({is: '('});
                    continue wordLoop;
                  }
                  parseState.unshift({
                    is: 'custom transform params',
                    transform: word,
                    on: 0,
                    params: [],
                    paramTypes: BUILT_IN_TRANSFORMS[word]
                  });
                  parseState.unshift({
                    is: parseState[0].paramTypes[0].type + ' param'
                  });
                  parseState.unshift({is: '('});
                  continue wordLoop;
                }
              }
              for(let transform in customFunctions) {
                if(word === transform){
                  console.log(`found ${transform}`);
                  if(customFunctions[word].length === 0){
                    parseState[0].transforms.push([word,[]]);
                    parseState.unshift({is: 'after transform'});
                    parseState.unshift({is: ')'});
                    parseState.unshift({is: '('});
                    continue wordLoop;
                  }
                  parseState.unshift({
                    is: 'custom transform params',
                    transform: word,
                    on: 0,
                    params: [],
                    paramTypes: customFunctions[word].params
                  });
                  parseState.unshift({
                    is: parseState[0].paramTypes[0].type + ' param'
                  });
                  parseState.unshift({is: '('});
                  continue wordLoop;
                }
              }
              newGeneralError(`Function undefined: ${word}`);
            break;
            default:
              newError('transform');
          }
        break;
        case 'after transform':
          switch(wordType){
            case '->':
              parseState.shift();
            break;
            case ';':
              switch(parseState[2].is){
                case 'choose items':
                  parseState[2].items.push([parseState[2].weight,parseState[1].transforms]);
                  delete parseState[2].weight;
                  parseState.shift();
                  parseState.shift();
                break;
                case 'transform':
                  parseState[2].transforms.push(parseState[1].transforms);
                  parseState.shift();
                  parseState.shift();
                break;
                case 'body':
                  parseState[3].body = parseState[1].transforms;
                  parseState.shift();
                  parseState.shift();
                  parseState.shift();
                break;
                case 'camera':
                  parseState[3].camera = parseState[1].transforms;
                  parseState.shift();
                  parseState.shift();
                  parseState.shift();
                break;
                default:
                  throw [`Unhandeled state transition from "${parseState[1].is}" to "${parseState[2].is}"`, parseState];
              }
            break;
            default:
            newError('";" or "->"');
          }
        break;
        case 'choose items':
          switch(wordType){
            case 'number':
              parseState[0].weight = parseFloat(word);
              parseState.unshift({is: 'transform', transforms: []});
              parseState.unshift({is: ':'});
            break;
            case '}':
              switch(parseState[1].is){
                case 'choose items':
                  parseState[1].items.push([parseState[1].weight,parseState[0].items]);
                  delete parseState[1].weight;
                  parseState.shift();
                break;
                case 'transform':
                  parseState[1].transforms.push(parseState[0].items);
                  parseState.shift();
                  parseState.unshift({is: 'after transform'});
                break;
                default:
                  throw [`Unhandeled state transition from "${parseState[0].is}" to "${parseState[1].is}"`, parseState];
              }
            break;
            default:
              newError('weight value or "}"');
          }
        break;
        default:
          throw [`Unhandeled state: ${parseState[0].is}`, parseState];
      }
    }
  }

  if(parseState.length === 1){
    delete parseState[0].is;

    if(!parseState[0].hasOwnProperty('body')){
      parseState[0].body = [];
    }

    if(!parseState[0].hasOwnProperty('camera')){
      parseState[0].camera = [];
    }

    parseState[0].customFunctions = customFunctions;

    errorModal.style.display = "none";
    //console.log('-- compiled --');
    return parseState[0];
  }
  General3arthError({at:code[code.length-2].length-1}, code.length-1, code[code.length-2])('unexpected end of input');
}

function from3arthLang(code) {
//console.log('-- compile --');
  try {
    return parseEverything(code);
  } catch (e) {
    console.error(e);
  }
}

let stuffToDo;

getFile("kleinian.3arth", r => {
  //console.log(r);
  stuffToDo = from3arthLang(r);
  //console.log(stuffToDo);
  refreshRender();
});
