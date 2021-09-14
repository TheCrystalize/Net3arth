function pointSymmetry(centerX, centerY, order) {
  return function(z) {
    let idr = Math.floor(Math.random() * order),
      dx = z.re - centerX, dy = z.im - centerY,
      da = (2 * Math.PI) / order,
      angle = idx * da,
      cosa = Math.cos(angle)
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
const pointSymmetryData = {
  name: "pointSymmetry",
  parameters: [{
      name: "centerX",
      type: "number",
      default: 0
    },
    {
      name: "centerY",
      type: "number",
      default: 0
    },
    {
      name: "order",
      type: "number",
      default: 1
    }
  ]
}
