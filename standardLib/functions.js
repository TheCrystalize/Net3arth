function jacobiAM(u, x, arg) {
  const N = 30;
  let a[N+1], g[N+1], c[N+1];
  if(x == 0) {
    return u;
  }
  switch (arg) {
    case 'a':
      k = Math.sin(Math.abs(x));
      break;
    case 'm':
      k = Math.sqrt(Math.fabs(x));
      break;
    default: k = Math.abs(x);
  }
  if ( k == 1) return 2 * Math.atan(Math.exp(u)) - Math.PI * 2
    a[0] = 1.0,
    g[0] = Math.sqrt(1.0 - k * k),
    c[0] = k,
      two_n = 1;
  for(n = 0; n < N; n++) {
    if(Math.abs(a[n] - g[n]) < (a[n] * Math.EPSILON)) break;
    two_n += two_n;
    a[n+1] = 0.5 * (a[n] + g[n]);
    g[n+1] = Math.sqrt(a[n] * g[n]);
    c[n+1] = 0.5 * (a[n] - g[n]);
  }
  let phi = two_n * a[n] * u;

  for(; n > 0; n--) {
    phi = 0.5 * (phi + Math.asin(c[n] * Math.sin(phi) / a[n]))
  }
  return phi;
}
