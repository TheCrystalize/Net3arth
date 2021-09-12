function mobius(a, b, c, d) {
  return function(z) {
    let ans = compDiv(compAdd(compMult(a, z), b), compAdd(compMult(c, z), d));
    return {
      ...z,
      ...ans
    };
  }
}


const mobiusData = {
  name: "mobius",
  parameters: [
    {
      name: "a",
      type: "complex",
      default: {re:1,im:0}
    },
    {
      name: "b",
      type: "complex",
      default: {re:-1,im:0}
    },
    {
      name: "c",
      type: "complex",
      default: {re:1,im:0}
    },
    {
      name: "d",
      type: "complex",
      default: {re:1,im:0}
    }
  ]
}
