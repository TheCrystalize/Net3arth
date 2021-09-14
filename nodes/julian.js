function julian(pow, dist) {
  return function(z) {
    let a = Math.atan2(z.im, z.re) + Math.floor(pow * Math.random()) * Math.PI * 2.0 / pow,
      r = Math.pow(z.re * z.re + z.im * z.im, dist / pow * 0.5);
    let ans = {
      re: Math.cos(a) * r,
      im: Math.sin(a) * r
    }
    return {
      ...z,
      ...ans
    }
  }
}
const julianData = {
  name: "julian",
  parameters: [{
      name: "pow",
      type: "float",
      default: 1
    },
    {
      name: "dist",
      type: "float",
      default: 1
    }
  ]
}
