/**
 * Busuanzi Visitor Counter - Modern ES6+ Version
 * A lightweight visitor statistics counter for websites
 *
 * Displays three types of statistics:
 * - site_pv: Total site page views
 * - page_pv: Current page views
 * - site_uv: Unique visitors to the site
 */

class BusuanziCounter {
  constructor(options = {}) {
    this.apiUrl = options.apiUrl || 'https://counter.busuanzi.icodeq.com/';
    this.statsTypes = ['site_pv', 'page_pv', 'site_uv'];
    this.retryAttempts = options.retryAttempts || 3;
    this.retryDelay = options.retryDelay || 1000;
  }

  /**
   * Initialize the counter when DOM is ready
   */
  async init() {
    if (document.readyState === 'loading') {
      document.addEventListener('DOMContentLoaded', () => this.fetchAndDisplay());
    } else {
      await this.fetchAndDisplay();
    }
  }

  /**
   * Fetch statistics using JSONP (required for cross-origin without CORS)
   */
  async fetchWithJsonp() {
    return new Promise((resolve, reject) => {
      const callbackName = `BusuanziCallback_${Math.floor(Math.random() * 1e16)}`;
      const timeout = setTimeout(() => {
        this.cleanup(callbackName, script);
        reject(new Error('Request timeout'));
      }, 10000);

      window[callbackName] = (data) => {
        clearTimeout(timeout);
        this.cleanup(callbackName, script);
        resolve(data);
      };

      const script = document.createElement('script');
      script.src = `${this.apiUrl}?jsonpCallback=${callbackName}`;
      script.referrerPolicy = 'no-referrer-when-downgrade';
      script.onerror = () => {
        clearTimeout(timeout);
        this.cleanup(callbackName, script);
        reject(new Error('Script loading failed'));
      };

      document.head.appendChild(script);
    });
  }

  /**
   * Fetch statistics using modern Fetch API (requires CORS support)
   */
  async fetchWithFetch() {
    try {
      const response = await fetch(this.apiUrl, {
        method: 'GET',
        mode: 'cors',
        credentials: 'omit',
        referrerPolicy: 'no-referrer-when-downgrade'
      });

      if (!response.ok) {
        throw new Error(`HTTP error! status: ${response.status}`);
      }

      return await response.json();
    } catch (error) {
      console.warn('Fetch API failed, falling back to JSONP:', error.message);
      return this.fetchWithJsonp();
    }
  }

  /**
   * Fetch and display statistics with retry logic
   */
  async fetchAndDisplay(attempt = 1) {
    try {
      // Try Fetch API first, with JSONP fallback
      const stats = await this.fetchWithFetch();
      this.updateDisplay(stats);
      this.showContainers();
    } catch (error) {
      console.error(`Busuanzi fetch attempt ${attempt} failed:`, error.message);

      if (attempt < this.retryAttempts) {
        await this.delay(this.retryDelay * attempt);
        await this.fetchAndDisplay(attempt + 1);
      } else {
        console.error('Busuanzi: All retry attempts failed');
        this.hideContainers();
      }
    }
  }

  /**
   * Update the display elements with fetched statistics
   */
  updateDisplay(stats) {
    this.statsTypes.forEach(type => {
      const element = document.getElementById(`busuanzi_value_${type}`);
      if (element && stats[type] !== undefined) {
        element.textContent = this.formatNumber(stats[type]);
      }
    });
  }

  /**
   * Format numbers with thousand separators
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
   * Hide all counter containers (on error)
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
   * Cleanup JSONP callback and script tag
   */
  cleanup(callbackName, script) {
    if (window[callbackName]) {
      delete window[callbackName];
    }
    if (script && script.parentElement) {
      script.parentElement.removeChild(script);
    }
  }

  /**
   * Delay helper for retry logic
   */
  delay(ms) {
    return new Promise(resolve => setTimeout(resolve, ms));
  }
}

// Auto-initialize with default options
const busuanzi = new BusuanziCounter();
busuanzi.init();

// Export for module usage
if (typeof module !== 'undefined' && module.exports) {
  module.exports = BusuanziCounter;
}
