---
# File Download Tips — v2
---

# File Download Tips

General utilities for downloading files from external services.

---

## Box Shared Links

Box shared links (e.g., `https://institution.app.box.com/s/<token>`) are **preview pages**, not direct download URLs.
A direct `curl` or `wget` to a `/s/<token>` URL will return a 404.

**Direct download URL pattern:**
```
https://institution.app.box.com/shared/static/<token>
```

Replace `/s/` with `/shared/static/` using the same token from the shared link.

**Command:**
```bash
curl -L -o output_filename.ext "https://institution.app.box.com/shared/static/<token>"
```

- `-L` is required to follow redirects (Box uses a redirect chain as part of its download flow).
- Works for large files — tested on a 15 GB `.Rds` object at ~100 MB/s.
