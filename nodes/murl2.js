function murl2(c, pow) {
  return function(z) {
    let angle = Math.atan2(z.im, z.re) * pow,
        cosa = Math.cos(angle), sina = Math.sin(angle),
        r = c * Math.pow(z.re * z.re + z.im * z.im, 0.5 * pow),
        real = r * cosa + 1, imag = r * sina,
        r1 = Math.pow(real * real + imag * imag, 1.0 / pow),
        angle1 = Math.atan2(imag, real) * 2 / pow,
        cosa1 = Math.cos(angle1), sina1 = Math.sin(angle1),
        re2 = r1 * cosa1, im2 = r1 * sina1,
        vp = c!= -1.0 ? Math.pow(c + 1, 2.0 / pow) : 0.0,
        r2 = vp / (r1 * r1);
    let ans = {
      re: (z.re * re2 + z.im * im2) * r2,
      im: (z.im * re2 - z.re * im2) * r2
    };
    return {
      ...z,
      ...ans
    };
  }
}


const murl2Data = {
  name: "murl2",
  parameters: [
    {
      name: "c",
      type: "number",
      default: 0
    },
    {
      name: "pow",
      type: "number",
      default: 2
    }
  ]
}
