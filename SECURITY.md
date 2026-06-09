# Security Policy

## Supported Versions

Only the latest release receives security fixes. Once a new version is published,
the previous release is no longer supported.

| Version | Supported |
| ------- | --------- |
| 1.0.x   | ✅ Yes    |
| < 1.0   | ❌ No     |

## Reporting a Vulnerability

**Please do not open a public GitHub issue for security vulnerabilities.**

Use GitHub's private vulnerability reporting feature instead:

1. Go to the [Security Advisories](https://github.com/wind-whispered/NumOptimizationkit/security/advisories) page.
2. Click **"Report a vulnerability"**.
3. Fill in the details described below.

If you are unable to use the above channel, you may contact the maintainer
directly via the email address on the GitHub profile
[@wind-whispered](https://github.com/wind-whispered).

### What to include in your report

- A clear description of the vulnerability and its potential impact.
- Wolfram Language / Mathematica version and operating system.
- NumOptimizationkit version (`PackageInformation["NumOptimizationkit"]`).
- A minimal, self-contained reproduction case (code snippet or notebook).
- Any suggested fix or mitigation, if you have one.

### Response timeline

| Milestone                              | Target     |
| -------------------------------------- | ---------- |
| Acknowledgement of report              | 3 days     |
| Initial assessment (valid / not valid) | 7 days     |
| Patch or mitigation published          | 30 days    |
| Public disclosure (coordinated)        | 90 days    |

We follow a **coordinated disclosure** policy: we ask reporters to allow us
reasonable time to investigate and release a fix before making details public.
We will credit you in the release notes unless you prefer to remain anonymous.

## Scope

NumOptimizationkit is a pure numerical-computing library with no network access,
no authentication logic, and no persistent storage. The realistic attack surface
is narrow, but the following classes of issues are considered **in scope**:

- **Arbitrary code execution** via unsafe evaluation of user-supplied symbolic
  expressions inside library code (e.g. unguarded `ToExpression` / `ReleaseHold`).
- **Unbounded resource consumption** (memory or CPU) triggered by crafted inputs
  that cause algorithms to diverge or allocate unexpectedly — when this behaviour
  is not documented and cannot be avoided by the caller.
- **Dependency vulnerabilities** in the paclet manifest or build toolchain that
  affect users who install the package.

The following are **out of scope**:

- Bugs that only affect results (numerical inaccuracy, wrong output, crashes on
  edge-case inputs) — please open a regular [GitHub issue](https://github.com/wind-whispered/NumOptimizationkit/issues) instead.
- Vulnerabilities in Wolfram Language / Mathematica itself — report these to
  [Wolfram Research](https://www.wolfram.com/legal/terms/wolfram-mathematica.html).
- Issues that require the attacker to already have arbitrary code execution on
  the target machine.

## License

This security policy applies to the source code distributed under the
[MIT License](LICENSE).
