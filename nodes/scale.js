function scale(s) {
  return function(z) {
    let ans = compMultScalar(z, s);
    return {
      ...z,
      ...ans
    }
  }
}

const scaleData = {
  name: "scale",
  parameters: [
    {
      name: "s",
      type: "float"
    }
  ]
}
