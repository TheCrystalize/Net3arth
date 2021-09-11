const mobius = {
  name: "mobius",
  function: function(a, b, c, d){
    return function(z){
      return compDiv(compAdd(compMult(a, z), c), compAdd(compMult(c, z), d));
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
