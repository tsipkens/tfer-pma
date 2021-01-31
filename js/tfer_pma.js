//==========================================================================//
// TFER_PMA CODE -----------------------------------------------------------//
var pi = 3.14159265359

var prop_pma = function(opts) {
  if (typeof opts == 'undefined' || (opts != null && opts.hasOwnProperty("__kwargtrans__"))) {
    ;
    var opts = 'Olfert';
  };
  if (opts == 'Olfert') {
    var prop = {
      'r1': 0.06,
      'r2': 0.061,
      'L': 0.2,
      'p': 1,
      'T': 293,
      'Q': (3 / 1000) / 60,
      'omega_hat': 32 / 33
    };
  } else if (opts == 'Buckley') {
    var prop = {
      'r1': 0.025,
      'r2': 0.024,
      'L': 0.1,
      'omega': ((13350 * 2) * pi) / 60,
      'p': 1,
      'T': 298,
      'Q': 0.00102 / 60,
      'omega_hat': 1
    };
  }
  prop['rc'] = (prop['r1'] + prop['r2']) / 2;
  prop['r_hat'] = prop['r1'] / prop['r2'];
  prop['del'] = (prop['r2'] - prop['r1']) / 2;
  prop['A'] = pi * (Math.pow(prop['r2'], 2) - Math.pow(prop['r1'], 2));
  prop['v_bar'] = prop['Q'] / prop['A'];
  var kB = 1.3806488e-23;
  prop['D'] = (function __lambda__(B) {
    return (kB * prop['T']) * B;
  });
  prop['Dm'] = 3;
  prop['m0'] = 4.7124e-25; // mass @ dm = 1 in fg
  return prop;
};


var get_setpoint = function(prop) {
  if (typeof prop == 'undefined' || (prop != null && prop.hasOwnProperty("__kwargtrans__"))) {
    ;
    var prop = {};
  };
  var sp = {
    'm_star': null,
    'V': null,
    'omega': null,
    'omega1': null,
    'Rm': null
  };
  if (arguments.length == 3) {
    sp[arguments[1]] = arguments[2];
    sp['Rm'] = 3;
  } else if (arguments.length == 5) {
    sp[arguments[1]] = arguments[2];
    sp[arguments[3]] = arguments[4];
  }
  var e = 1.60218e-19;
  if (sp['m_star'] == null) {
    if (!(sp['omega1'])) {
      sp['omega1'] = sp['omega'] / ((Math.pow(prop['r_hat'], 2) - prop['omega_hat']) /
        (Math.pow(prop['r_hat'], 2) - 1) + ((Math.pow(prop['r1'], 2) * (prop['omega_hat'] - 1)) /
        (Math.pow(prop['r_hat'], 2) - 1)) / Math.pow(prop['rc'], 2));
    }
    sp['alpha'] = (sp['omega1'] * (Math.pow(prop['r_hat'], 2) - prop['omega_hat'])) /
      (Math.pow(prop['r_hat'], 2) - 1);
    sp['beta'] = ((sp['omega1'] * Math.pow(prop['r1'], 2)) * (prop['omega_hat'] - 1)) /
      (Math.pow(prop['r_hat'], 2) - 1);
    sp['m_star'] = sp['V'] / ((Math.log(1 / prop['r_hat']) / e) *
      Math.pow(sp['alpha'] * prop['rc'] + sp['beta'] / prop['rc'], 2));
    sp['omega'] = sp['alpha'] + sp['beta'] / Math.pow(prop['rc'], 2);
    sp['omega2'] = sp['alpha'] + sp['beta'] / Math.pow(prop['r2'], 2);
  } else if (sp['omega1'] != null) {
    sp['alpha'] = (sp['omega1'] * (Math.pow(prop['r_hat'], 2) - prop['omega_hat'])) /
      (Math.pow(prop['r_hat'], 2) - 1);
    sp['beta'] = ((sp['omega1'] * Math.pow(prop['r1'], 2)) * (prop['omega_hat'] - 1)) /
      (Math.pow(prop['r_hat'], 2) - 1);
    sp['V'] = ((sp['m_star'] * Math.log(1 / prop['r_hat'])) / e) *
      Math.pow(sp['alpha'] * prop['rc'] + sp['beta'] / prop['rc'], 2);
    sp['omega2'] = sp['alpha'] + sp['beta'] / Math.pow(prop['r2'], 2);
    sp['omega'] = sp['alpha'] + sp['beta'] / Math.pow(prop['rc'], 2);
  } else if (sp['omega'] != null) {
    sp['omega1'] = sp['omega'] / ((Math.pow(prop['r_hat'], 2) - prop['omega_hat']) /
      (Math.pow(prop['r_hat'], 2) - 1) + ((Math.pow(prop['r1'], 2) * (prop['omega_hat'] - 1)) /
      (Math.pow(prop['r_hat'], 2) - 1)) / Math.pow(prop['rc'], 2));
    sp['alpha'] = (sp['omega1'] * (Math.pow(prop['r_hat'], 2) - prop['omega_hat'])) /
      (Math.pow(prop['r_hat'], 2) - 1);
    sp['beta'] = ((sp['omega1'] * Math.pow(prop['r1'], 2)) * (prop['omega_hat'] - 1)) /
      (Math.pow(prop['r_hat'], 2) - 1);
    sp['V'] = ((sp['m_star'] * Math.log(1 / prop['r_hat'])) / e) *
      Math.pow(sp['alpha'] * prop['rc'] + sp['beta'] / prop['rc'], 2);
    sp['omega2'] = sp['alpha'] + sp['beta'] / Math.pow(prop['r2'], 2);
  } else if (sp['V'] != null) {
    var v_theta_rc = Math.sqrt((sp['V'] * e) / (sp['m_star'] * Math.log(1 / prop['r_hat'])));
    var A = (prop['rc'] * (Math.pow(prop['r_hat'], 2) - prop['omega_hat'])) /
      (Math.pow(prop['r_hat'], 2) - 1) + (1 / prop['rc']) * ((Math.pow(prop['r1'], 2) *
        (prop['omega_hat'] - 1)) / (Math.pow(prop['r_hat'], 2) - 1));
    sp['omega1'] = v_theta_rc / A;
    sp['alpha'] = (sp['omega1'] * (Math.pow(prop['r_hat'], 2) - prop['omega_hat'])) /
      (Math.pow(prop['r_hat'], 2) - 1);
    sp['beta'] = ((sp['omega1'] * Math.pow(prop['r1'], 2)) * (prop['omega_hat'] - 1)) /
      (Math.pow(prop['r_hat'], 2) - 1);
    sp['omega2'] = sp['alpha'] + sp['beta'] / Math.pow(prop['r2'], 2);
    sp['omega'] = sp['alpha'] + sp['beta'] / Math.pow(prop['rc'], 2);
  } else if (sp['Rm'] != null) {
    var n_B = get_nb(sp['m_star'], prop);
    var __left0__ = mp2zp(sp['m_star'], 1, prop['T'], prop['p'], prop);
    var B_star = __left0__[0];
    sp['m_max'] = sp['m_star'] * (1 / sp['Rm'] + 1);
    sp['omega'] = Math.sqrt(prop['Q'] / ((((((sp['m_star'] * B_star) * 2) * pi) *
      Math.pow(prop['rc'], 2)) * prop['L']) * (Math.pow(sp['m_max'] / sp['m_star'], n_B + 1) -
      Math.pow(sp['m_max'] / sp['m_star'], n_B))));
    sp['omega1'] = sp['omega'] / ((Math.pow(prop['r_hat'], 2) - prop['omega_hat']) /
      (Math.pow(prop['r_hat'], 2) - 1) + ((Math.pow(prop['r1'], 2) * (prop['omega_hat'] - 1)) /
      (Math.pow(prop['r_hat'], 2) - 1)) / Math.pow(prop['rc'], 2));
    sp['alpha'] = (sp['omega1'] * (Math.pow(prop['r_hat'], 2) - prop['omega_hat'])) /
      (Math.pow(prop['r_hat'], 2) - 1);
    sp['beta'] = ((sp['omega1'] * Math.pow(prop['r1'], 2)) * (prop['omega_hat'] - 1)) /
      (Math.pow(prop['r_hat'], 2) - 1);
    sp['omega2'] = sp['alpha'] + sp['beta'] / Math.pow(prop['r2'], 2);
    sp['V'] = ((sp['m_star'] * Math.log(1 / prop['r_hat'])) / e) * Math.pow(sp['alpha'] * prop['rc'] +
      sp['beta'] / prop['rc'], 2);
  } else {
    console.log('No setpoint specified.');
  }
  var m_star = sp['m_star'];
  if (!(sp['Rm'])) {
    var __left0__ = get_resolution(sp['m_star'], sp['omega'], prop);
    sp['Rm'] = __left0__[0];
    sp['m_max'] = __left0__[1];
  }
  return [sp, m_star];
};

var get_resolution = function(m_star, omega, prop) {
  print('Finding resolution...');
  var n_B = get_nb(m_star, prop);
  var __left0__ = mp2zp(m_star, 1, prop['T'], prop['p'], prop);
  var B_star = __left0__[0];
  var t0 = prop['Q'] / ((((((m_star * B_star) * 2) * pi) * prop['L']) * Math.pow(omega, 2)) * Math.pow(prop['rc'], 2));
  var m_rat = (function __lambda__(Rm) {
    return 1 / Rm + 1;
  });
  var fun = (function __lambda__(Rm) {
    return Math.pow(m_rat(Rm), n_B + 1) - Math.pow(m_rat(Rm), n_B);
  });
  var Rm = scipy.optimize.fmin((function __lambda__(Rm) {
    return Math.pow(t0 - fun(Rm), 2);
  }), __kwargtrans__({
    x0: 5
  }));
  var Rm = Rm[0];
  var m_max = m_star * (1 / Rm + 1);
  print('Complete.');
  print(' ');
  return [Rm, m_max];
};

var get_nb = function(m_star, prop) {
  var m_high = m_star * 1.001;
  var m_low = m_star * 0.999;
  var __left0__ = mp2zp(m_high, 1, prop['T'], prop['p'], prop);
  var B_high = __left0__[0];
  var __left0__ = mp2zp(m_low, 1, prop['T'], prop['p'], prop);
  var B_low = __left0__[0];
  var n_B = Math.log10(B_high / B_low) / Math.log10(m_high / m_low);
  return n_B;
};

var dm2zp = function(d, z, T, p) {
  if (typeof T == 'undefined' || (T != null && T.hasOwnProperty("__kwargtrans__"))) {
    ;
    var T = 0.0;
  };
  if (typeof p == 'undefined' || (p != null && p.hasOwnProperty("__kwargtrans__"))) {
    ;
    var p = 0.0;
  };
  var e = 1.6022e-19;
  if (T == 0.0 || p == 0.0) {
    var mu = 1.82e-05;
    var B = Cc(d) / (((3 * pi) * mu) * d);
  } else {
    var S = 110.4;
    var T_0 = 296.15;
    var vis_23 = 1.83245e-05;
    var mu = (vis_23 * Math.pow(T / T_0, 1.5)) * ((T_0 + S) / (T + S));
    var B = Cc(d, T, p) / (((3 * pi) * mu) * d);
  }
  var Zp = (B * e) * z;
  return [B, Zp];
};

var Cc = function(d, T, p) {
  if (typeof T == 'undefined' || (T != null && T.hasOwnProperty("__kwargtrans__"))) {
    ;
    var T = 0.0;
  };
  if (typeof p == 'undefined' || (p != null && p.hasOwnProperty("__kwargtrans__"))) {
    ;
    var p = 0.0;
  };
  if (T == 0.0 || p == 0.0) {
    var mfp = 6.65e-08;
    var A1 = 1.257;
    var A2 = 0.4;
    var A3 = 0.55;
  } else {
    var S = 110.4;
    var mfp_0 = 6.73e-08;
    var T_0 = 296.15;
    var p_0 = 101325;
    var p = p * p_0;
    var mfp = ((mfp_0 * Math.pow(T / T_0, 2)) * (p_0 / p)) * ((T_0 + S) / (T + S));
    var A1 = 1.165;
    var A2 = 0.483;
    var A3 = 0.997 / 2;
  }
  var Kn = (2 * mfp) / d;
  var Cc = 1 + Kn * (A1 + A2 * Math.exp(-(2 * A3) / Kn));
  return Cc;
};

var mp2zp = function(m, z, T, P, prop) {
  if (typeof T == 'undefined' || (T != null && T.hasOwnProperty("__kwargtrans__"))) {
    ;
    var T = 0.0;
  };
  if (typeof P == 'undefined' || (P != null && P.hasOwnProperty("__kwargtrans__"))) {
    ;
    var P = 0.0;
  };
  if (typeof prop == 'undefined' || (prop != null && prop.hasOwnProperty("__kwargtrans__"))) {
    ;
    var prop = {};
  };
  if (!('m0' in prop) || !('Dm' in prop)) {
    sys.exit('Please specify the mass-mobility relation parameters in prop.');
  }
  var d = Math.pow(m / prop['m0'], 1 / prop['Dm']) * 1e-9;
  if (T == 0.0 || P == 0.0) {
    var __left0__ = dm2zp(d, z);
    var B = __left0__[0];
    var Zp = __left0__[1];
  } else {
    var __left0__ = dm2zp(d, z, T, P);
    var B = __left0__[0];
    var Zp = __left0__[1];
  }
  return [B, Zp, d];
};


var tfer_1S = function(sp, m, d, z, prop) {
  if (typeof prop == 'undefined' || (prop != null && prop.hasOwnProperty("__kwargtrans__"))) {
    ;
    var prop = {};
  };
  var __left0__ = parse_inputs(sp, m, d, z, prop);
  var tau = __left0__[0];
  var rs = __left0__[3];
  var lam = (((2 * tau) * (Math.pow(sp['alpha'], 2) - Math.pow(sp['beta'], 2) / Math.pow(rs, 4))) * prop['L']) / prop['v_bar'];
  var G0 = (function __lambda__(r) {
    return rs + (r - rs) * Math.exp(-(lam));
  });
  var ra = Math.min(prop['r2'], Math.max(prop['r1'], G0(prop['r1'])));
  var rb = Math.min(prop['r2'], Math.max(prop['r1'], G0(prop['r2'])));
  var Lambda = (1 / (2 * prop['del'])) * (rb - ra);
  return [Lambda, G0];
};

var tfer_1C = function(sp, m, d, z, prop) {
  if (typeof prop == 'undefined' || (prop != null && prop.hasOwnProperty("__kwargtrans__"))) {
    ;
    var prop = {};
  };
  var __left0__ = parse_inputs(sp, m, d, z, prop);
  var tau = __left0__[0];
  var C0 = __left0__[1];
  var C3 = tau * (((Math.pow(sp['alpha'], 2) * prop['rc'] + ((2 * sp['alpha']) * sp['beta']) / prop['rc']) + Math.pow(sp['beta'], 2) / Math.pow(prop['rc'], 3)) - C0 / (m * prop['rc']));
  var C4 = tau * (((Math.pow(sp['alpha'], 2) - ((2 * sp['alpha']) * sp['beta']) / Math.pow(prop['rc'], 2)) - (3 * Math.pow(sp['beta'], 2)) / Math.pow(prop['rc'], 4)) + C0 / (m * Math.pow(prop['rc'], 2)));
  var G0 = (function __lambda__(r) {
    return (prop['rc'] + ((r - prop['rc']) + C3 / C4) * Math.exp((-(C4) * prop['L']) / prop['v_bar'])) - C3 / C4;
  });
  var ra = Math.min(prop['r2'], Math.max(prop['r1'], G0(prop['r1'])));
  var rb = Math.min(prop['r2'], Math.max(prop['r1'], G0(prop['r2'])));
  var Lambda = (1 / (2 * prop['del'])) * (rb - ra);
  return [Lambda, G0];
};

var tfer_1C_diff = function(sp, m, d, z, prop) {
  if (typeof prop == 'undefined' || (prop != null && prop.hasOwnProperty("__kwargtrans__"))) {
    ;
    var prop = {};
  };
  var __left0__ = parse_inputs(sp, m, d, z, prop);
  var D = __left0__[2];
  var sig = Math.sqrt(((2 * prop['L']) * D) / prop['v_bar']);
  var __left0__ = tfer_1C(sp, m, d, z, prop);
  var G0 = __left0__[1];
  var rho_fun = (function __lambda__(G, r) {
    return (G - r) / (Math.sqrt(2) * sig);
  });
  var kap_fun = (function __lambda__(G, r) {
    return (G - r) * math.erf(rho_fun(G, r)) + (sig * Math.sqrt(2 / pi)) * Math.exp(-(Math.pow(rho_fun(G, r), 2)));
  });
  var K22 = kap_fun(G0(prop['r2']), prop['r2']);
  var K21 = kap_fun(G0(prop['r2']), prop['r1']);
  var K12 = kap_fun(G0(prop['r1']), prop['r2']);
  var K11 = kap_fun(G0(prop['r1']), prop['r1']);
  var Lambda = (-(1) / (4 * prop['del'])) * (((K22 - K12) - K21) + K11);
  Lambda[K22 > 100.0] = 0;
  Lambda[Math.abs(Lambda) < 1e-10] = 0;
  return [Lambda, G0];
};

var tfer_W1 = function(sp, m, d, z, prop) {
  var __left0__ = parse_inputs(sp, m, d, z, prop);
  var tau = __left0__[0];
  var C0 = __left0__[1];
  var rs = __left0__[3];
  var lam = (((2 * tau) * (Math.pow(sp['alpha'], 2) - Math.pow(sp['beta'], 2) / Math.pow(rs, 4))) * prop['L']) / prop['v_bar'];
  var G0 = (function __lambda__(r) {
    return (1 / (sp['omega1'] * Math.sqrt(m))) * Math.sqrt(((m * Math.pow(sp['omega1'], 2)) * Math.pow(r, 2) - C0) * Math.exp(-(lam)) + C0);
  });
  var ra = Math.min(prop['r2'], Math.max(prop['r1'], G0(prop['r1'])));
  var rb = Math.min(prop['r2'], Math.max(prop['r1'], G0(prop['r2'])));
  var Lambda = (1 / (2 * prop['del'])) * (rb - ra);
  return [Lambda, G0];
};

var tfer_W1_diff = function(sp, m, d, z, prop) {
  var __left0__ = parse_inputs(sp, m, d, z, prop);
  var D = __left0__[2];
  var sig = Math.sqrt(((2 * prop['L']) * D) / prop['v_bar']);
  var __left0__ = tfer_W1(sp, m, d, z, prop);
  var G0 = __left0__[1];
  var rho_fun = (function __lambda__(G, r) {
    return (G - r) / (Math.sqrt(2) * sig);
  });
  var kap_fun = (function __lambda__(G, r) {
    return (G - r) * math.erf(rho_fun(G, r)) + (sig * Math.sqrt(2 / pi)) * Math.exp(-(Math.pow(rho_fun(G, r), 2)));
  });
  var K22 = kap_fun(G0(prop['r2']), prop['r2']);
  var K21 = kap_fun(G0(prop['r2']), prop['r1']);
  var K12 = kap_fun(G0(prop['r1']), prop['r2']);
  var K11 = kap_fun(G0(prop['r1']), prop['r1']);
  var Lambda = (-(1) / (4 * prop['del'])) * (((K22 - K12) - K21) + K11);
  Lambda[K22 > 100.0] = 0;
  Lambda[Math.abs(Lambda) < 1e-10] = 0;
  return [Lambda, G0];
};

var parse_inputs = function(sp, m, d, z, prop) {
  if (typeof d == 'undefined' || (d != null && d.hasOwnProperty("__kwargtrans__"))) {
    ;
    var d = 0;
  };
  if (typeof z == 'undefined' || (z != null && z.hasOwnProperty("__kwargtrans__"))) {
    ;
    var z = 1;
  };
  if (typeof prop == 'undefined' || (prop != null && prop.hasOwnProperty("__kwargtrans__"))) {
    ;
    var prop = {};
  };
  var e = 1.60218e-19;
  var q = z * e;
  if (d == 0 || d == null) {
    console.log('Invoking mass-mobility relation to determine Zp.');
    var __left0__ = mp2zp(m, z, prop['T'], prop['p'], prop);
    var B = __left0__[0];
  } else {
    var __left0__ = dm2zp(d, z, prop['T'], prop['p']);
    var B = __left0__[0];
  }
  var tau = B * m;
  var D = prop['D'](B) * z;
  var C0 = (sp['V'] * q) / Math.log(1 / prop['r_hat']);

  // Note: Whether to pick the +ive of -ive root for rs is chosen based on a
  // heuristic approach. Specifically, the root closer to the centerline is
  // chosen, except when the -ive root is zero (which is the case for APM
  // conditions, where the +ive root should always be used).
  // evaluate +ive and -ive roots
  r_m = (math.sqrt(C0 / m) - math.sqrt(C0 / m - 4 * sp['alpha'] * sp['beta'])) / (2 * sp['alpha']);
  r_p = (math.sqrt(C0 / m) + math.sqrt(C0 / m - 4 * sp['alpha'] * sp['beta'])) / (2 * sp['alpha']);

  // determine which root is closer to centerline radius
  function indexOfSmallest(a) {
    var lowest = 0;
    for (var i = 1; i < a.length; i++) {
      if (a[i] < a[lowest]) lowest = i;
    }
    return lowest;
  }
  idx = indexOfSmallest([Math.abs(r_m - prop['rc']), Math.abs(r_p - prop['rc'])]);
  if (r_m == 0) {
    idx = 2; // avoid zero values for APM case
  }

  // assign one of the two roots to rs
  rs = r_m; // by default use -ive root
  if (idx == 2) {
    rs = r_p // if closer to +ive root, use +ive root
  }

  // zero out cases where no equilibrium radius exists (also removes complex numbers)
  if (C0 / m < (4 * sp['alpha'] * sp['beta'])) {
    rs = 0;
  }

  rs = Math.max(rs, 1e-20) // avoids division by zero for z = 0 case

  return [tau, C0, D, rs];
};

var linspace = function(a, b, n) {
  if (typeof n === "undefined") n = Math.max(Math.round(b - a) + 1, 1);
  if (n < 2) {
    return n === 1 ? [a] : [];
  }
  var i, ret = Array(n);
  n--;
  for (i = n; i >= 0; i--) {
    ret[i] = (i * b + (n - i) * a) / n;
  }
  return ret;
};


//-- CHARGING FUNCTIONS ----------------------------------------//
//   Incomplete at the moment.
var tfer_charge = function(d, z) {

  var e = 1.602177e-19, // elementary charge
    epi = 8.85418e-12, // dielectric constant (for air) [F/m]
    kB = 1.38065e-23, // Boltzmann's constant
    Z_Z = 0.875, // ion mobility ratio (Wiedensohler, 1988)
    T = 298; // temperature

  var fn = Array(z.length).fill(0).map(x => Array(d.length).fill(0)); // initialize fn

  for (dd in d) {
    for (zz in z) {
      if (z[zz] < 3) { // if charge state less than 3
        // Gopalakrishnan below z = 3
        a = [
          [-0.3880, -8.0157, -40.714],
          [0.4545, 3.2536, 17.487],
          [-0.1634, -0.5018, -2.6146],
          [0.0091, 0.0223, 0.1282]
        ];
        exponent = 0;
        for (jj = 0; jj < 4; jj++) { // loop through coefficients in a
          exponent = exponent + (a[jj][z[zz]] * Math.log(d[dd] * 1e9) ** jj);
        }
        fn[zz][dd] = Math.exp(exponent);
      } else { // Wiedensohler for z = 3 and above
        fn[zz][dd] = e / Math.sqrt(4 * pi * pi * epi * kB * T * d[dd]) *
          Math.exp(0 - (z[zz] - 2 * pi * epi * kB * T * Math.log(Z_Z) * d[dd] / e ** 2) ** 2 /
            (4 * pi * epi * kB * T * d[dd] / e ** 2));

        if (fn[zz][dd] < 6e-5) {
          fn[zz][dd] = 0
        } // truncate small values
      }
    }
  }

  return fn
}


var fCharge = 0;
var parse_fun = function(sp, m, d, prop, fun) {
  var Lambda = Array(m.length),
      tCharge = tfer_charge(d, z_vec); // unused as y-scale becomes challenging
  for (ii in m) { // loop over particle mass
    Lambda[ii] = 0 // initialize at zero
    for (zz in z_vec) { // loop over integer charge states
      __left0__ = fun(sp, m[ii], d[ii], z_vec[zz], prop)
      if (fCharge==1) {
         Lambda[ii] = Lambda[ii] + __left0__[0] * tCharge[zz][ii]
      } else {
        Lambda[ii] = Lambda[ii] + __left0__[0]
      }
    }
  }
  return Lambda;
};


console.log('Load complete.')
console.log(' ')

// END TFER_PMA CODE -------------------------------------------------------//
//==========================================================================//
