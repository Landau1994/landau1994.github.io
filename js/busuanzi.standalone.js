/**
 * Standalone Visitor Counter - No External API Required
 * Uses localStorage for tracking and optional custom backend support
 *
 * Features:
 * - Works completely offline
 * - Privacy-friendly (all data stays local)
 * - Can sync with custom backend if needed
 * - Tracks page views and unique visitors
 */

class StandaloneCounter {
  constructor(options = {}) {
    this.storageKey = options.storageKey || 'visitor_stats';
    this.cookieKey = options.cookieKey || 'visitor_id';
    this.backend = options.backend || null; // Optional backend URL
    this.statsTypes = ['site_pv', 'page_pv', 'site_uv'];
    this.sessionKey = 'session_pages_viewed';
  }

  /**
   * Initialize the counter when DOM is ready
   */
  async init() {
    if (document.readyState === 'loading') {
      document.addEventListener('DOMContentLoaded', () => this.run());
    } else {
      await this.run();
    }
  }

  /**
   * Main execution flow
   */
  async run() {
    try {
      // Get or create visitor ID
      const visitorId = this.getOrCreateVisitorId();

      // Get current page URL
      const pageUrl = this.getCurrentPage();

      // Load stats from localStorage
      let stats = this.loadStats();

      // Update stats
      stats = this.updateStats(stats, visitorId, pageUrl);

      // Save stats
      this.saveStats(stats);

      // Display stats
      this.displayStats(stats, pageUrl);

      // Optionally sync with backend
      if (this.backend) {
        await this.syncWithBackend(stats);
      }

      this.showContainers();
    } catch (error) {
      console.error('StandaloneCounter error:', error);
      this.hideContainers();
    }
  }

  /**
   * Get or create unique visitor ID using cookies
   */
  getOrCreateVisitorId() {
    // Check if visitor ID exists in cookie
    let visitorId = this.getCookie(this.cookieKey);

    if (!visitorId) {
      // Generate new visitor ID
      visitorId = this.generateId();
      // Set cookie for 1 year
      this.setCookie(this.cookieKey, visitorId, 365);
    }

    return visitorId;
  }

  /**
   * Generate unique ID
   */
  generateId() {
    return `${Date.now()}-${Math.random().toString(36).substr(2, 9)}`;
  }

  /**
   * Get current page identifier (pathname)
   */
  getCurrentPage() {
    return window.location.pathname || '/';
  }

  /**
   * Load stats from localStorage
   */
  loadStats() {
    try {
      const data = localStorage.getItem(this.storageKey);
      if (data) {
        return JSON.parse(data);
      }
    } catch (error) {
      console.warn('Failed to load stats from localStorage:', error);
    }

    // Default stats structure
    return {
      site_pv: 0,      // Total page views across site
      site_uv: 0,      // Unique visitors
      visitors: {},    // Visitor tracking
      pages: {}        // Per-page statistics
    };
  }

  /**
   * Update statistics
   */
  updateStats(stats, visitorId, pageUrl) {
    // Increment total site page views
    stats.site_pv++;

    // Track unique visitor
    if (!stats.visitors[visitorId]) {
      stats.visitors[visitorId] = {
        firstVisit: Date.now(),
        visitCount: 0,
        pages: []
      };
      stats.site_uv++;
    }

    // Update visitor info
    stats.visitors[visitorId].visitCount++;
    stats.visitors[visitorId].lastVisit = Date.now();

    // Track page in visitor's history
    if (!stats.visitors[visitorId].pages.includes(pageUrl)) {
      stats.visitors[visitorId].pages.push(pageUrl);
    }

    // Track per-page statistics
    if (!stats.pages[pageUrl]) {
      stats.pages[pageUrl] = {
        pv: 0,
        uv: 0,
        visitors: []
      };
    }

    stats.pages[pageUrl].pv++;

    // Track unique visitor for this page
    if (!stats.pages[pageUrl].visitors.includes(visitorId)) {
      stats.pages[pageUrl].visitors.push(visitorId);
      stats.pages[pageUrl].uv++;
    }

    return stats;
  }

  /**
   * Save stats to localStorage
   */
  saveStats(stats) {
    try {
      localStorage.setItem(this.storageKey, JSON.stringify(stats));
    } catch (error) {
      console.error('Failed to save stats:', error);
    }
  }

  /**
   * Display statistics on page
   */
  displayStats(stats, pageUrl) {
    // Site-wide page views
    this.updateElement('site_pv', stats.site_pv);

    // Site-wide unique visitors
    this.updateElement('site_uv', stats.site_uv);

    // Current page views
    const pageStats = stats.pages[pageUrl] || { pv: 0 };
    this.updateElement('page_pv', pageStats.pv);
  }

  /**
   * Update a display element
   */
  updateElement(type, value) {
    const element = document.getElementById(`busuanzi_value_${type}`);
    if (element) {
      element.textContent = this.formatNumber(value);
    }
  }

  /**
   * Format number with thousand separators
   */
  formatNumber(num) {
    return num.toString().replace(/\B(?=(\d{3})+(?!\d))/g, ',');
  }

  /**
   * Show all counter containers
   */
  showContainers() {
    this.statsTypes.forEach(type => {
      const container = document.getElementById(`busuanzi_container_${type}`);
      if (container) {
        container.style.display = 'inline';
      }
    });
  }

  /**
   * Hide all counter containers
   */
  hideContainers() {
    this.statsTypes.forEach(type => {
      const container = document.getElementById(`busuanzi_container_${type}`);
      if (container) {
        container.style.display = 'none';
      }
    });
  }

  /**
   * Cookie management
   */
  setCookie(name, value, days) {
    const expires = new Date();
    expires.setTime(expires.getTime() + days * 24 * 60 * 60 * 1000);
    document.cookie = `${name}=${value};expires=${expires.toUTCString()};path=/;SameSite=Lax`;
  }

  getCookie(name) {
    const nameEQ = name + "=";
    const ca = document.cookie.split(';');
    for (let i = 0; i < ca.length; i++) {
      let c = ca[i];
      while (c.charAt(0) === ' ') c = c.substring(1, c.length);
      if (c.indexOf(nameEQ) === 0) return c.substring(nameEQ.length, c.length);
    }
    return null;
  }

  /**
   * Optional: Sync with custom backend
   */
  async syncWithBackend(stats) {
    if (!this.backend) return;

    try {
      const response = await fetch(this.backend, {
        method: 'POST',
        headers: {
          'Content-Type': 'application/json'
        },
        body: JSON.stringify({
          action: 'track',
          page: this.getCurrentPage(),
          timestamp: Date.now()
        })
      });

      if (response.ok) {
        const serverStats = await response.json();
        // Optionally update with server stats
        if (serverStats.site_pv) {
          this.updateElement('site_pv', serverStats.site_pv);
        }
        if (serverStats.site_uv) {
          this.updateElement('site_uv', serverStats.site_uv);
        }
        if (serverStats.page_pv) {
          this.updateElement('page_pv', serverStats.page_pv);
        }
      }
    } catch (error) {
      console.warn('Backend sync failed:', error.message);
      // Silently fail - local stats still work
    }
  }

  /**
   * Get current statistics (for debugging)
   */
  getStats() {
    return this.loadStats();
  }

  /**
   * Reset all statistics (for testing)
   */
  reset() {
    localStorage.removeItem(this.storageKey);
    console.log('Statistics reset');
  }
}

// Auto-initialize
const standaloneCounter = new StandaloneCounter();
standaloneCounter.init();

// Export for module usage and debugging
if (typeof module !== 'undefined' && module.exports) {
  module.exports = StandaloneCounter;
}

// Expose to window for debugging
if (typeof window !== 'undefined') {
  window.StandaloneCounter = StandaloneCounter;
  window.counterInstance = standaloneCounter;
}
