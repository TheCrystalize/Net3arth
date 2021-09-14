function splits(x, y) {
  return function(z) {
    let xoff = z.re >= 0 ? x : -x,
      yoff = z.im >= 0 ? y : -y;
    let ans = {
      re: z.re + xoff,
      im: z.im + yoff
    }
    return {
      ...z,
      ...ans
    }
  }
}

const splitsData = {
  name: "splits",
  parameters: [
    {
      name: "x",
      type: "number",
      default: 0
    },
    {
      name: "y",
      type: "number",
      default: 0
    }
  ]
}
