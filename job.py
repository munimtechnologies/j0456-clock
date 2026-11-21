# Copyright 2026 Munim, Inc.
# SPDX-License-Identifier: Apache-2.0
"""Run the pre-registered J0456-20 test on public Swift archive data.

    python job.py --data-dir /data --out-dir /out

1. Lists every Swift observation within 0.2 deg of eRASSt J045650.3-203750 (HEASARC TAP).
2. Downloads the cleaned XRT photon-counting event files and UVOT sky images.
3. Measures 0.3-2 keV aperture counts per observation; checks the cycles Liu et al. 2024
   already published (validation gates); then grades the 2023-2026 data exactly as
   PREREGISTRATION.md says.
Made with AI agents working for MV2 of Munim, Inc.
"""
import argparse
import csv
import json
import os
import re
import statistics
import sys
import time
import urllib.error
import urllib.parse
import urllib.request

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

from jz import __version__, ephem, geometry, grade, photometry, uvot  # noqa: E402

# Host nucleus, WISEA J045649.80-203748.0 (the eROSITA/XMM position agrees to a few arcsec).
RA = (4 + 56 / 60.0 + 49.80 / 3600.0) * 15.0
DEC = -(20 + 37 / 60.0 + 48.0 / 3600.0)
CONE_DEG = 0.2
MAX_OFFAXIS_ARCMIN = 8.0
TAP = "https://heasarc.gsfc.nasa.gov/xamin/vo/tap/sync"
ARCHIVE = "https://heasarc.gsfc.nasa.gov/FTP/swift/data/obs"
UA = "python-urllib/3 research-client"


def log(*a):
    print(*a, flush=True)


def http_get(url, tries=5):
    last = None
    for i in range(tries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": UA})
            with urllib.request.urlopen(req, timeout=180) as r:
                return r.read()
        except urllib.error.HTTPError as e:
            if e.code == 404:
                return None
            last = e
        except Exception as e:  # network hiccup: retry
            last = e
        time.sleep(3 * (i + 1))
    raise RuntimeError("GET failed after %d tries: %s (%s)" % (tries, url, last))


def download(url, dest):
    if os.path.exists(dest) and os.path.getsize(dest) > 0:
        return dest
    data = http_get(url)
    if data is None:
        return None
    os.makedirs(os.path.dirname(dest), exist_ok=True)
    tmp = dest + ".part"
    with open(tmp, "wb") as f:
        f.write(data)
    os.replace(tmp, dest)
    time.sleep(0.2)
    return dest


def list_dir(url):
    data = http_get(url)
    if data is None:
        return None
    return re.findall(r'href="([^"?/]+)"', data.decode("latin-1"))


def query_observations(out_dir):
    q = ("SELECT obsid,start_time,xrt_exposure,uvot_exposure,ra,dec FROM swiftmastr "
         "WHERE CONTAINS(POINT('ICRS',ra,dec),CIRCLE('ICRS',%.6f,%.6f,%.3f))=1 "
         "ORDER BY start_time" % (RA, DEC, CONE_DEG))
    url = TAP + "?" + urllib.parse.urlencode(
        {"REQUEST": "doQuery", "LANG": "ADQL", "FORMAT": "text", "QUERY": q})
    raw = http_get(url)
    if raw is None:
        raise RuntimeError("HEASARC TAP query returned 404")
    text = raw.decode("utf-8", "replace")
    with open(os.path.join(out_dir, "swiftmastr_query.txt"), "w", encoding="utf-8") as f:
        f.write(q + "\n\n" + text)
    obs = []
    for line in text.splitlines():
        parts = [p.strip() for p in line.split("|")]
        if len(parts) >= 6 and re.fullmatch(r"\d{11}", parts[0]):
            def num(s):
                try:
                    return float(s)
                except ValueError:
                    return 0.0
            obs.append({"obsid": parts[0], "mjd": num(parts[1]), "xrt_exposure": num(parts[2]),
                        "uvot_exposure": num(parts[3])})
    if not obs:
        raise RuntimeError("no observations parsed from the TAP response")
    obs.sort(key=lambda o: o["mjd"])
    return obs


def find_obs_dir(obsid, mjd, sub):
    for shift in (0, -1, 1):
        ym = geometry.mjd_to_archive_month(mjd, shift)
        base = "%s/%s/%s/" % (ARCHIVE, ym, obsid)
        names = list_dir(base + sub)
        if names is not None:
            return base + sub, names
    return None, None


def read_pc_events(path):
    from astropy.io import fits
    from astropy.wcs import WCS

    with fits.open(path, memmap=False) as h:
        ev = h["EVENTS"]
        hdr = ev.header
        names = [c.upper() for c in ev.columns.names]
        ix = names.index("X") + 1
        iy = names.index("Y") + 1
        if not str(hdr["TCTYP%d" % ix]).startswith("RA---TAN"):
            raise RuntimeError("unexpected X axis type in %s" % path)
        w = dict(crval1=float(hdr["TCRVL%d" % ix]), crval2=float(hdr["TCRVL%d" % iy]),
                 crpix1=float(hdr["TCRPX%d" % ix]), crpix2=float(hdr["TCRPX%d" % iy]),
                 cdelt1=float(hdr["TCDLT%d" % ix]), cdelt2=float(hdr["TCDLT%d" % iy]))
        d = ev.data
        if d is None or len(d) == 0:
            xs, ys, pis = [], [], []
        else:
            xs = [float(v) for v in d["X"]]
            ys = [float(v) for v in d["Y"]]
            pis = [int(v) for v in d["PI"]]
        exposure = float(hdr.get("EXPOSURE") or hdr.get("LIVETIME") or hdr.get("ONTIME") or 0.0)
        ra_pnt = hdr.get("RA_PNT")
        dec_pnt = hdr.get("DEC_PNT")
    x0, y0 = geometry.tan_world_to_pixel(RA, DEC, **w)
    aw = WCS(naxis=2)
    aw.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    aw.wcs.crval = [w["crval1"], w["crval2"]]
    aw.wcs.crpix = [w["crpix1"], w["crpix2"]]
    aw.wcs.cdelt = [w["cdelt1"], w["cdelt2"]]
    ax, ay = aw.all_world2pix([[RA, DEC]], 1)[0]
    if abs(ax - x0) > 0.5 or abs(ay - y0) > 0.5:
        raise RuntimeError("projection mismatch %s: ours (%.2f,%.2f) astropy (%.2f,%.2f)"
                           % (path, x0, y0, ax, ay))
    pixscale = abs(w["cdelt1"]) * 3600.0
    return xs, ys, pis, x0, y0, pixscale, exposure, ra_pnt, dec_pnt


def xrt_for_obs(o, data_dir):
    url, names = find_obs_dir(o["obsid"], o["mjd"], "xrt/event/")
    if url is None:
        return {"status": "NO_DIR"}
    pat = re.compile(r"sw%sxpcw\dpo_cl\.evt\.gz" % o["obsid"])
    pc = sorted(n for n in names if pat.fullmatch(n))
    if not pc:
        return {"status": "NO_PC_FILE"}
    n_src = n_bkg = 0
    exposure = 0.0
    offaxis = []
    for n in pc:
        p = download(url + n, os.path.join(data_dir, "xrt", o["obsid"], n))
        if p is None:
            continue
        xs, ys, pis, x0, y0, pixscale, exp, ra_pnt, dec_pnt = read_pc_events(p)
        if exp <= 0.0:
            continue
        s, b = photometry.count_aperture(xs, ys, pis, x0, y0, pixscale)
        n_src += s
        n_bkg += b
        exposure += exp
        if ra_pnt is not None and dec_pnt is not None:
            offaxis.append(geometry.angsep_arcsec(RA, DEC, float(ra_pnt), float(dec_pnt)) / 60.0)
    if exposure <= 0.0:
        return {"status": "NO_EXPOSURE"}
    off = max(offaxis) if offaxis else None
    status = "OK"
    if off is not None and off > MAX_OFFAXIS_ARCMIN:
        status = "OFFAXIS"
    return {"status": status, "n_src": n_src, "n_bkg": n_bkg, "exposure_s": exposure,
            "offaxis_arcmin": off, "files": ";".join(pc)}


def uvot_for_obs(o, data_dir):
    url, names = find_obs_dir(o["obsid"], o["mjd"], "uvot/image/")
    if url is None:
        return []
    pat = re.compile(r"sw%su(w1|w2|m2)_sk\.img\.gz" % o["obsid"])
    out = []
    for n in sorted(names):
        m = pat.fullmatch(n)
        if not m:
            continue
        p = download(url + n, os.path.join(data_dir, "uvot", o["obsid"], n))
        if p is None:
            continue
        try:
            r = uvot.phot_sky_image(p, RA, DEC)
        except Exception as e:  # a bad image must not stop the X-ray test
            log("  uvot error", n, e)
            continue
        if r is None:
            continue
        r.update({"obsid": o["obsid"], "mjd": o["mjd"], "filter": m.group(1),
                  "ab_mag": uvot.ab_mag(r["net_rate"], m.group(1))})
        out.append(r)
    return out


def uv_secondary(uv_rows):
    """Descriptive only: median UVOT rate per filter in the published and forecast windows."""
    wins = {"cycle4_plateau": [grade.C4_PLATEAU], "cycle4_faint": [grade.C4_FAINT],
            "cycle6_on": [grade.ON6], "cycle5_6_off": [grade.OFF5, grade.OFF6]}
    res = {}
    for f in ("m2", "w2", "w1"):
        fr = [r for r in uv_rows if r["filter"] == f]
        item = {}
        for name, ws in wins.items():
            sel = [r["net_rate"] for r in fr if any(w[0] <= r["mjd"] <= w[1] for w in ws)]
            item[name] = {"n": len(sel),
                          "median_rate": statistics.median(sel) if sel else None}
        res[f] = item
    return res


def write_csv(path, rows, cols):
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=cols, extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow(r)


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-dir", default="data")
    ap.add_argument("--out-dir", default="out")
    ap.add_argument("--skip-uvot", action="store_true")
    args = ap.parse_args(argv)
    os.makedirs(args.out_dir, exist_ok=True)
    os.makedirs(args.data_dir, exist_ok=True)
    t0 = time.time()

    obs = query_observations(args.out_dir)
    log("observations in cone:", len(obs))

    xrt_rows = []
    skipped = []
    cols = ["obsid", "mjd", "status", "exposure_s", "n_src", "n_bkg", "rate", "rate_err",
            "sigma", "ul3_rate", "detected", "level", "forecast_phase", "offaxis_arcmin", "files"]
    for i, o in enumerate(obs):
        if o["xrt_exposure"] <= 0.0:
            skipped.append({"obsid": o["obsid"], "mjd": o["mjd"], "status": "NO_XRT_EXPOSURE"})
            continue
        try:
            r = xrt_for_obs(o, args.data_dir)
        except Exception as e:
            r = {"status": "ERROR: %s" % e}
        r.update({"obsid": o["obsid"], "mjd": o["mjd"]})
        if r["status"] == "OK":
            r.update(photometry.measure(r["n_src"], r["n_bkg"], r["exposure_s"]))
            r["forecast_phase"] = ephem.phase_class(o["mjd"])
            xrt_rows.append(r)
        else:
            skipped.append(r)
        log("%3d/%d %s mjd=%.2f %s src=%s bkg=%s exp=%s" % (
            i + 1, len(obs), o["obsid"], o["mjd"], r["status"], r.get("n_src"), r.get("n_bkg"),
            r.get("exposure_s")))
        if (i + 1) % 20 == 0:
            write_csv(os.path.join(args.out_dir, "xrt_obs.csv"), xrt_rows, cols)
    write_csv(os.path.join(args.out_dir, "xrt_obs.csv"), xrt_rows, cols)
    write_csv(os.path.join(args.out_dir, "xrt_skipped.csv"), skipped,
              ["obsid", "mjd", "status", "offaxis_arcmin", "exposure_s", "n_src", "n_bkg"])

    gates = grade.run_gates(xrt_rows)
    result = grade.grade(xrt_rows)
    if not gates["all_passed"]:
        result["outcome_if_gates_had_passed"] = result["outcome"]
        result["outcome"] = "UNGRADED_GATE_FAILURE"

    uv_rows = []
    if not args.skip_uvot:
        for o in obs:
            if o["uvot_exposure"] > 0.0:
                try:
                    uv_rows.extend(uvot_for_obs(o, args.data_dir))
                except Exception as e:
                    log("  uvot failure", o["obsid"], e)
        write_csv(os.path.join(args.out_dir, "uvot_obs.csv"), uv_rows,
                  ["obsid", "mjd", "filter", "net_rate", "net_rate_err", "ab_mag", "exposure_s",
                   "n_ext", "coi_factor", "src_counts", "bkg_counts_in_src"])
    uv = uv_secondary(uv_rows) if uv_rows else {"status": "NOT_RUN"}

    out = {
        "code_version": __version__,
        "target": {"name": "eRASSt J045650.3-203750", "ra_deg": RA, "dec_deg": DEC},
        "forecast": {"source": "Liu et al. 2024, A&A 683, L13 (arXiv:2401.14091)",
                     "recurrence_days": ephem.RECUR, "rise_start_cycle5": ephem.T5,
                     "predicted_rise6": ephem.rise_start(6),
                     "predicted_drop6": ephem.predicted_drop(6)},
        "method": {"band_kev": [0.3, 2.0], "src_radius_arcsec": photometry.SRC_R,
                   "bkg_annulus_arcsec": [photometry.BKG_IN, photometry.BKG_OUT],
                   "psf_fraction": photometry.EEF, "max_offaxis_arcmin": MAX_OFFAXIS_ARCMIN},
        "n_observations_in_cone": len(obs),
        "n_xrt_measured": len(xrt_rows),
        "n_skipped": len(skipped),
        "gates": gates,
        "primary": result,
        "secondary_uv": uv,
        "runtime_s": round(time.time() - t0, 1),
    }
    with open(os.path.join(args.out_dir, "results.json"), "w", encoding="utf-8") as f:
        json.dump(out, f, indent=1, default=str)

    lines = ["# eRASSt J0456-20 after May 2023: results (code %s)" % __version__, "",
             "Gates: %s (all passed: %s)" % (json.dumps(gates["status"]), gates["all_passed"]),
             "", "PRIMARY OUTCOME: **%s**" % result["outcome"],
             "z(on vs off) = %s" % result.get("z_on_vs_off"),
             "depth reading: %s" % result.get("depth_reading"),
             "T2 amplitude: %s" % json.dumps(result["T2_amplitude"], default=str),
             "T3 timing: %s" % result["T3_timing"]["status"]]
    for k in ("on_stack", "off_stack", "post_stack"):
        s = result.get(k)
        if s:
            lines.append("%s: n_obs=%d exp=%.0f s src=%d bkg=%d rate=%.4g +/- %.2g sigma=%.2f "
                         "ul3=%.4g detected=%s" % (k, s["n_obs"], s["exposure_s"], s["n_src"],
                                                   s["n_bkg"], s["rate"], s["rate_err"],
                                                   s["sigma"], s["ul3_rate"], s["detected"]))
    for g in result["T4_late_epochs"]:
        lines.append("late group %.1f-%.1f: n=%d exp=%.0f rate=%.4g +/- %.2g ul3=%.4g det=%s "
                     "phases=%s" % (g["mjd_first"], g["mjd_last"], g["n_obs"], g["exposure_s"],
                                    g["rate"], g["rate_err"], g["ul3_rate"], g["detected"],
                                    ",".join(g["forecast_phases"])))
    lines.append("")
    lines.append("Secondary UV: %s" % json.dumps(uv, default=str))
    with open(os.path.join(args.out_dir, "summary.md"), "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")
    log("\n".join(lines))
    return 0


if __name__ == "__main__":
    sys.exit(main())
