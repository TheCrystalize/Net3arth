function parseVal(txt) {
  let v = 0;
  try {
    v = eval(txt);
    if(typeof v !== 'number') {
      throw 'NaN';
    }
  } catch (e) {
    v = parseFloat(txt);
  }
  return isNaN(v) ? 0 : v;
}

function loadUI() {
  
}
