function translate(c){
  return function(z){
    let ans = compAdd(z, c);
    return {
      ...z,
      ...ans
    };
  }
}

const translateData = {
  name: "translate",
  parameters: [
    {
      name: "c",
      type: "complex",
      default: {
        re: 0,
        im: 0
      }
    }
  ]
}
