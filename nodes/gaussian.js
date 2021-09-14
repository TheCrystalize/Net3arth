function gaussian(pow) {
  let a = Math.random() * 2 * Math.PI,
      r = gaussRnd() * pow;
  return function(z) {
    let ans = {
      re: Math.cos(a) * r + z.re,
      im: Math.sin(a) * r + z.im
    };
    return {
      ...z,
      ...ans
    };
  }
}


const gaussianData = {
  name: "gaussian",
  parameters: [
    {
      name: "pow",
      type: "number",
      default: 1.0
    }
  ]
}
