# Busuanzi Counter Versions - Comparison Guide

## 📦 Available Versions

### 1. **busuanzi.pure.mini.js** (Original)
- ❌ Relies on external API (counter.busuanzi.icodeq.com)
- ❌ Currently experiencing downtime
- ✅ Minified for production
- ✅ Cross-site statistics

### 2. **busuanzi.modern.js** (Modern API Version)
- ❌ Still relies on external API
- ✅ Clean ES6+ code
- ✅ Fetch API with JSONP fallback
- ✅ Better error handling
- ✅ Retry logic

### 3. **busuanzi.standalone.js** (Recommended)
- ✅ **No external API required**
- ✅ Works completely offline
- ✅ Privacy-friendly (all data local)
- ✅ Instant loading (no network delay)
- ✅ 100% reliable (no API downtime)
- ⚠️ Statistics are per-device (not shared)
- ✅ Optional backend sync available

---

## 🚀 Quick Start - Standalone Version

### Basic Usage (No Backend)

```html
<!-- Include the script -->
<script src="js/busuanzi.standalone.js"></script>

<!-- Add counter elements anywhere -->
<span id="busuanzi_container_site_pv">
  Total views: <span id="busuanzi_value_site_pv"></span>
</span>

<span id="busuanzi_container_page_pv">
  Page views: <span id="busuanzi_value_page_pv"></span>
</span>

<span id="busuanzi_container_site_uv">
  Visitors: <span id="busuanzi_value_site_uv"></span>
</span>
```

### With Optional Backend

```html
<script src="js/busuanzi.standalone.js"></script>
<script>
  // Override with custom backend
  const counter = new StandaloneCounter({
    backend: 'https://your-api.com/track'
  });
  counter.init();
</script>
```

---

## 🧪 Testing

### Test the Standalone Version

Open in your browser:
```
http://localhost:8888/busuanzi-standalone-test.html
```

**What you'll see:**
- ✅ Instant counter display (no API wait)
- ✅ Numbers increment on refresh
- ✅ Works in incognito (new visitor)
- ✅ Debug tools and controls

---

## 🔧 Optional Backend Setup

If you want to sync statistics across devices/users:

1. **Install dependencies:**
```bash
pip install flask flask-cors
```

2. **Run the backend:**
```bash
python3 busuanzi-backend.py
```

3. **Update your frontend:**
```javascript
const counter = new StandaloneCounter({
  backend: 'http://localhost:5000/api/track'
});
counter.init();
```

---

## 📊 Feature Comparison

| Feature | Original | Modern | **Standalone** |
|---------|----------|--------|----------------|
| External API | Required | Required | **Optional** |
| Works Offline | ❌ | ❌ | **✅** |
| Privacy | Low | Low | **High** |
| Reliability | Low | Medium | **High** |
| Loading Speed | Slow | Slow | **Instant** |
| Cross-device | ✅ | ✅ | With backend |
| ES6+ Code | ❌ | ✅ | **✅** |
| Error Handling | ❌ | ✅ | **✅** |

---

## 🎯 Which Version to Use?

### Use **Standalone** (Recommended) if:
- ✅ You want reliable, instant counters
- ✅ Privacy is important
- ✅ You don't need cross-device sync
- ✅ You have a static site (GitHub Pages, Netlify, etc.)

### Use **Standalone + Backend** if:
- ✅ You need cross-device statistics
- ✅ You have server infrastructure
- ✅ You want centralized analytics

### Use **Modern** if:
- ⚠️ You must use Busuanzi API (when it's back online)
- ⚠️ You need cross-site statistics

---

## 🔨 Migration from Original to Standalone

**Good news:** It's a drop-in replacement!

1. Replace the script:
```html
<!-- OLD -->
<script src="js/busuanzi.pure.mini.js"></script>

<!-- NEW -->
<script src="js/busuanzi.standalone.js"></script>
```

2. HTML stays the same! No changes needed to:
- `busuanzi_container_*` divs
- `busuanzi_value_*` spans

---

## 🐛 Debugging

```javascript
// View current statistics
window.counterInstance.getStats()

// Reset statistics (for testing)
window.counterInstance.reset()

// Check visitor ID
document.cookie
```

---

## 📝 Notes

- **localStorage** stores up to 5-10MB (plenty for statistics)
- **Cookies** track unique visitors (expires in 1 year)
- **No personal data** is collected or transmitted
- **Works in all modern browsers** (Chrome, Firefox, Safari, Edge)

---

## 🎉 Summary

The **standalone version** is the best choice for most use cases:
- No external dependencies
- Privacy-friendly
- Fast and reliable
- Easy to deploy

Test it now: http://localhost:8888/busuanzi-standalone-test.html
