# Two-Phase Microflow Simulator

A browser-based phase-field solver for immiscible laminar two-phase flow (water and cyclohexane by default), set up with COMSOL's *Two-Phase Flow, Phase Field* formulation so the two can be compared.

## What's in this folder

| File | What it is | Do you edit it? |
|---|---|---|
| `index.html` | The website: layout, inputs, plots and the methods write-up. It loads `solver.js`. | Yes, for text, layout, presets and styling |
| `solver.js` | The physics engine: Navier–Stokes and Cahn–Hilliard on a staggered grid. Used by both the website and the tests. | Yes, for the physics and numerics |
| `assets/lab-background.jpg` | Background and banner image | Swap it for another image with the same name |
| `tests/test_poiseuille.js` | Validation: channel flow against the exact parabola | Run it after changing the solver |
| `tests/test_droplets.js` | Validation: Laplace pressure and sessile-drop contact angle | Run it after changing the solver |
| `.nojekyll` | Empty file that tells GitHub Pages to serve the files as they are | No |

## Your edit → check → deploy loop

1. **Edit** `solver.js` or `index.html` in VS Code and save.
2. **Look at it locally.** Right-click `index.html` and choose *Open with Live Server* (install the free "Live Server" extension once). The page reloads by itself every time you save.
3. **Check the physics** after a solver change. In VS Code's terminal (Ctrl + `):
   ```
   node tests/test_poiseuille.js
   node tests/test_droplets.js laplace
   node tests/test_droplets.js sessile 60
   ```
   The errors should stay close to the published numbers: Poiseuille about 0.25 % / 0.5 %, Laplace about 2 %, sessile about 58° for a 60° target.
4. **Deploy.** In VS Code's *Source Control* panel (the branch icon), type a short message such as "tweak mobility", click **Commit**, then **Sync / Push**. GitHub Pages republishes in about a minute.
5. If the live site still shows the old version, hard-refresh with Ctrl + Shift + R (Cmd + Shift + R on a Mac).

## One-time setup (GitHub Pages)

1. Create a free account at github.com. Make a **New repository** (for example `microflow-simulator`, Public).
2. On the new repository's page, click **uploading an existing file**. Drag in everything from this folder, **including** the `assets` and `tests` folders and `.nojekyll`, then commit.
3. Open **Settings → Pages**. Under *Build and deployment* choose **Deploy from a branch**, then branch `main` and folder `/ (root)`. Save.
4. After about a minute the site is live at `https://<your-username>.github.io/microflow-simulator/`.
5. In VS Code: press **Ctrl + Shift + P**, run **Git: Clone**, paste the repository URL, and open the cloned folder. That folder is the one you edit from now on, and steps 1–5 of the loop above deploy it.

Tip: on the repository page on github.com, pressing the `.` key opens a full VS Code editor in the browser, which is handy on a computer without your setup.
