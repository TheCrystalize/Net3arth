function mobius(a, b, c, d) {
  return function(z) {
    let ans = div(add(mult(a, z), b), add(mult(c, z), d));
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
      default: {re:0,im:0}
    },
    {
      name: "c",
      type: "complex",
      default: {re:0,im:0}
    },
    {
      name: "d",
      type: "complex",
      default: {re:1,im:0}
    }
  ]
}
