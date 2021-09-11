const mobius = {
  name: "mobius",
  function: function(a, b, c, d){
    return function(z){
      let ans = compDiv(compAdd(compMult(a, z), b), compAdd(compMult(c, z), d));
      return {
        ...z,
        ...ans
      }
    }
  },
  parameters: [
    {
      name: "a",
      type: "complex"
    },
    {
      name: "b",
      type: "complex"
    },
    {
      name: "c",
      type: "complex"
    },
    {
      name: "d",
      type: "complex"
    }
  ]
}
