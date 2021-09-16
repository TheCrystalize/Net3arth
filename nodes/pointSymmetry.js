function pointSymmetry(centerX, centerY, order) {
  return function(z) {
    let idr = Math.floor(Math.random() * order),
      dx = z.re - centerX, dy = z.im - centerY,
      da = (2 * Math.PI) / order;
    let angle = idr * da;
    let cosa = Math.cos(angle),
        sina = Math.sin(angle)
    }
    return {
      ...z,
      re: centerX + dx * cosa + dy * sina,
      im: centerY + dy * cosa - dx * sina
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
