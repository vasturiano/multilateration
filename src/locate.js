import { minimize_Powell } from 'optimization-js';
import minimize_cobyla from './3rdparty/jscobyla/index';

import { Earth, Point, Circle } from './geometry';

const sum = vals => vals.reduce((agg, v) => agg + v, 0);

const Norm = (x, y, mode='2d') => {
  if (mode === '2d') return ((x[0] - y[0]) ** 2 + (x[1] - y[1]) ** 2) ** .5;
  if (mode === '3d') return ((x[0] - y[0]) ** 2 + (x[1] - y[1]) ** 2 + (x[2] - y[2]) ** 2) ** .5;
  if (mode === 'earth') return Earth.gcd(x[0], x[1], y[0], y[1]);
};

const sumError = (x, c, r, mode) => {
  const l = c.length;
  let e = 0;
  for (let i = 0; i < l; i++) {
    e += (Norm(x, c[i].std(), mode) - r[i]) ** 2;
  }
  return e;
};

const isDisjoint = (cA, fg = 0) => {
  // does not support sophisticated checking for disjoint area on earth surface models
  const l = cA.length;
  for (let i = 0; i < l; i++) {
    for (let j = i+1; j < l; j++) {
      if (!cA[j].touch(cA[i], fg)) return true;
    }
  }
  return false;
};

function lse(cA, mode='2d', cons = true) {
  const l = cA.length;
  const r = cA.map(w => w.r);
  const c = cA.map(w => w.c);
  const S = sum(r);
  const W = r.map(w => (S - w) / ((l - 1) * S));

  let p0 = new Point(0, 0, 0); // Initialized point
  for (let i = 0; i < l; i++) {
    p0 = p0.sum(W[i]).sum(c[i]);
  }

  const x0 = (mode === '2d' || mode === 'earth')
    ? [p0.x, p0.y]
    : [p0.x, p0.y.p0.z];

  const fg1 = mode === 'earth' ? 1 : 0;

  let ans;
  if (cons) {
    // console.info('GC-LSE geolocating...');
    if (!isDisjoint(cA, fg1)) {
      const constrain = (x, beaconIdx) => r[beaconIdx] - Norm(x, c[beaconIdx].std(), mode);

      let x = x0.slice(); // result

      const res = minimize_cobyla(
        (n, m, x, con) => { // objective function
          for (let i = 0; i < m; i++) {
            con[i] = constrain(x, i); // one constrain per beacon
          }
          return sumError(x, c, r, mode);
        },
        x.length, // # vars
        l, // # constraints
        x0.slice(), // result
        1, // rhobeg
        1e-5, // rhoend
        0, // iprint
        1000 // maxfun
      );

      ans = x;
    } else {
      throw new Error('Disjoint');
    }
  } else {
    // console.info('LSE geolocating...');

    const res = minimize_Powell(x => sumError(x, c, r, mode), x0);

    ans = res.argument;
  }

  return new Point(...ans);
}

function locate(beacons, { mode = '2d', constrain = true }) {
  if (mode !== '2d' && mode !== '3d' && mode !== 'earth') {
    throw new Error(`Mode not supported: ${mode}`);
  }

  const circles = beacons.map(({ distance, ...coords }) => {
    const c = mode === '2d'
      ? new Point(coords.x, coords.y)
      : mode === '3d'
        ? new Point(coords.x, coords.y, coords.z)
        : new Point(coords.lng, coords.lat);

    return new Circle(c, distance); // in Earth mode, gcd distances are specified in meters
  });

  const { x, y, z } = lse(circles, mode, constrain);

  if (mode === 'earth') {
    let lng = x;
    let lat = y;

    // make lat, lngs uniform
    while (lng > 180) { lng -= 360 }
    while (lng < -180) { lng += 360 }
    while (lat > 90) { lat -= 180 }
    while (lat < -90) { lat += 180 }

    return { lat, lng };
  }

  return mode === '2d'
    ? { x, y }
    : { x, y, z };
}

export default locate;
