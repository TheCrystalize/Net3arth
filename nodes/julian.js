const julian = {
  name: "julian",
  function: function(pow, dist){
    return function(z){
      let a = Math.atan2(z.im, z.re) + Math.floor(pow * Math.random()) * Math.PI * 2.0 / pow
      let r = Math.pow(z.re * z.re + z.im * z.im, dist / pow * 0.5)
      let ans = {
        re: Math.cos(a) * r,
        im: Math.sin(a) * r
      }
      return {
        ...z,
        ...ans
      }
    }
  },
  parameters: [
    {
      name: "pow",
      type: "float"
    },
    {
      name: "dist",
      type: "float"
    }
  ]
}
