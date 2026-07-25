# Signing + notarizing the macOS FASTag desktop app

The self-contained `.app` is **built and verified** (see `gui/scripts/bundle-macos.sh`):
a single bundle carrying the CLI, its full dylib closure with **one** `libomp`
(which fixes `OMP: Error #15`), `share/OpenMS`, the ~1.1 GB taxonomy, and the
icon. It runs species detection standalone with no `KMP_DUPLICATE_LIB_OK`.

What's left is **Developer ID signing + notarization**, which needs Apple
credentials that must live in CI, never on a dev machine. This is the runbook.

## Secrets (add in GitHub → Settings → Secrets and variables → Actions)

Same set BALL/BALLView uses (see the `software-signing` runbook):

| secret | what |
|---|---|
| `MACOS_CERTIFICATE_BASE64` | Developer ID Application cert + key, `.p12`, base64 |
| `MACOS_CERTIFICATE_PASSWORD` | the `.p12` password |
| `MACOS_KEYCHAIN_PASSWORD` | password for the throwaway CI keychain |
| `MACOS_SIGNING_IDENTITY` | `Developer ID Application: Name (TEAMID)` |
| `MACOS_APPLE_ID` | Apple ID for notarytool |
| `MACOS_TEAM_ID` | Apple Developer Team ID |
| `MACOS_NOTARY_PASSWORD` | app-specific password for notarization |

## CI steps (add to `ci.yml`'s macOS leg, AFTER the CLI `dist/FASTag` is built,
## gated so it skips when the cert secret is absent — never breaks a release)

1. **Assemble the app's resources with the same closure the CLI build produced:**
   ```bash
   FASTAG_BIN=build/FASTag \
   OPENMS_LIB="$OPENMS_INSTALL/lib" DEPS_LIB="$CONDA_PREFIX/lib" \
   OPENMS_SHARE="$OPENMS_INSTALL/share/OpenMS" \
   TAXONOMY_DIR="$RUNNER_TEMP/taxonomy" \
   gui/scripts/bundle-macos.sh
   ```
   The taxonomy k-mer index isn't in the repo; fetch it once from the release
   (`gh release download "$TAG" -p 'FASTag-taxonomy-k7.tar.gz'`) into
   `$RUNNER_TEMP/taxonomy`, or rebuild it with `buildtaxdb`.

2. **Import the cert into a throwaway keychain** (decode `MACOS_CERTIFICATE_BASE64`
   → `.p12` → `security import` → set as default). Tear it down after.

3. **Sign every nested Mach-O inside out** (the CLI + all `lib/*.dylib`), hardened
   runtime + secure timestamp — notarization rejects any unsigned executable:
   ```bash
   find gui/src-tauri/resources/fastag \( -name '*.dylib' -o -name FASTag \) -print0 \
     | xargs -0 -I{} codesign --force --timestamp --options runtime \
         --sign "$MACOS_SIGNING_IDENTITY" {}
   ```

4. **Build + sign + notarize the app in one shot** — Tauri does the outer signing
   and notarization from env vars:
   ```bash
   cd gui && npm ci
   APPLE_CERTIFICATE="$MACOS_CERTIFICATE_BASE64" \
   APPLE_CERTIFICATE_PASSWORD="$MACOS_CERTIFICATE_PASSWORD" \
   APPLE_SIGNING_IDENTITY="$MACOS_SIGNING_IDENTITY" \
   APPLE_ID="$MACOS_APPLE_ID" APPLE_PASSWORD="$MACOS_NOTARY_PASSWORD" \
   APPLE_TEAM_ID="$MACOS_TEAM_ID" \
   npm run tauri build   # emits a signed, notarized, stapled .dmg + .app
   ```
   `tauri.conf.json` already sets `bundle.macOS.entitlements` →
   `entitlements.plist` (disables library validation for the bundled third-party
   dylibs) and `minimumSystemVersion`.

5. **Upload** `gui/src-tauri/target/release/bundle/dmg/FASTag_*.dmg` to the
   release with `gh release upload`.

## Gotchas (these bite on the first run)

- **`secrets.*` is illegal in `if:`** — map the cert secret into a step output
  and gate on that, or the whole workflow fails with a 0-job `startup_failure`.
- **Notarization is iterate-to-green.** The first submissions usually fail on a
  missed nested binary; `xcrun notarytool log <id>` names the exact file. "Sign
  every Mach-O", step 3, is why — don't hand-list suspects.
- **Two arches.** Build `macos-14` (arm64) and `macos-13` (x64) legs; the `.app`
  is ~1.5 GB each because the taxonomy is inside, so notarization upload is slow.
- Ad-hoc signatures from `bundle-macos.sh` (`codesign --sign -`) are placeholders
  the Developer ID pass in step 3/4 overwrites.
