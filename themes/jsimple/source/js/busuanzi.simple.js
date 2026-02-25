/**
 * Busuanzi 纯前端计数器 - 无需服务器
 * 完全在浏览器中运行，使用 localStorage 存储统计数据
 *
 * 特点：
 * - 无需任何后端服务器
 * - 数据存储在浏览器本地
 * - 即时加载，无网络延迟
 * - 隐私友好，数据不外传
 */

class BusuanziCounter {
  constructor() {
    this.storageKey = 'busuanzi_stats';
    this.cookieKey = 'busuanzi_visitor_id';
    this.statsTypes = ['site_pv', 'page_pv', 'site_uv'];
  }

  // 初始化计数器
  init() {
    if (document.readyState === 'loading') {
      document.addEventListener('DOMContentLoaded', () => this.run());
    } else {
      this.run();
    }
  }

  // 主执行流程
  run() {
    try {
      const visitorId = this.getOrCreateVisitorId();
      const pageUrl = this.getCurrentPage();
      let stats = this.loadStats();

      stats = this.updateStats(stats, visitorId, pageUrl);
      this.saveStats(stats);
      this.displayStats(stats, pageUrl);
      this.showContainers();
    } catch (error) {
      console.error('Busuanzi 计数器错误:', error);
      this.hideContainers();
    }
  }

  // 获取或创建访客ID
  getOrCreateVisitorId() {
    let visitorId = this.getCookie(this.cookieKey);
    if (!visitorId) {
      visitorId = `v_${Date.now()}_${Math.random().toString(36).substr(2, 9)}`;
      this.setCookie(this.cookieKey, visitorId, 365);
    }
    return visitorId;
  }

  // 获取当前页面路径
  getCurrentPage() {
    return window.location.pathname || '/';
  }

  // 从 localStorage 加载统计数据
  loadStats() {
    try {
      const data = localStorage.getItem(this.storageKey);
      if (data) {
        return JSON.parse(data);
      }
    } catch (error) {
      console.warn('无法从 localStorage 加载统计:', error);
    }
    return {
      site_pv: 0,
      site_uv: 0,
      visitors: {},
      pages: {}
    };
  }

  // 更新统计数据
  updateStats(stats, visitorId, pageUrl) {
    // 增加全站浏览量
    stats.site_pv++;

    // 跟踪唯一访客
    if (!stats.visitors[visitorId]) {
      stats.visitors[visitorId] = {
        firstVisit: Date.now(),
        visitCount: 0
      };
      stats.site_uv++;
    }
    stats.visitors[visitorId].visitCount++;
    stats.visitors[visitorId].lastVisit = Date.now();

    // 跟踪页面统计
    if (!stats.pages[pageUrl]) {
      stats.pages[pageUrl] = {
        pv: 0,
        uv: 0,
        visitors: []
      };
    }
    stats.pages[pageUrl].pv++;

    if (!stats.pages[pageUrl].visitors.includes(visitorId)) {
      stats.pages[pageUrl].visitors.push(visitorId);
      stats.pages[pageUrl].uv++;
    }

    return stats;
  }

  // 保存统计数据到 localStorage
  saveStats(stats) {
    try {
      localStorage.setItem(this.storageKey, JSON.stringify(stats));
    } catch (error) {
      console.error('无法保存统计数据:', error);
    }
  }

  // 显示统计数据
  displayStats(stats, pageUrl) {
    this.updateElement('site_pv', stats.site_pv);
    this.updateElement('site_uv', stats.site_uv);

    const pageStats = stats.pages[pageUrl] || { pv: 0 };
    this.updateElement('page_pv', pageStats.pv);
  }

  // 更新页面元素
  updateElement(type, value) {
    const element = document.getElementById(`busuanzi_value_${type}`);
    if (element) {
      element.textContent = this.formatNumber(value);
    }
  }

  // 格式化数字（添加千位分隔符）
  formatNumber(num) {
    return num.toString().replace(/\B(?=(\d{3})+(?!\d))/g, ',');
  }

  // 显示所有计数器容器
  showContainers() {
    this.statsTypes.forEach(type => {
      const container = document.getElementById(`busuanzi_container_${type}`);
      if (container) {
        container.style.display = 'inline';
      }
    });
  }

  // 隐藏所有计数器容器
  hideContainers() {
    this.statsTypes.forEach(type => {
      const container = document.getElementById(`busuanzi_container_${type}`);
      if (container) {
        container.style.display = 'none';
      }
    });
  }

  // Cookie 操作
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
      while (c.charAt(0) === ' ') c = c.substring(1);
      if (c.indexOf(nameEQ) === 0) return c.substring(nameEQ.length);
    }
    return null;
  }

  // 调试方法：获取当前统计
  getStats() {
    return this.loadStats();
  }

  // 调试方法：重置统计
  reset() {
    localStorage.removeItem(this.storageKey);
    console.log('统计数据已重置');
  }
}

// 自动初始化
const busuanzi = new BusuanziCounter();
busuanzi.init();

// 导出供调试使用
if (typeof window !== 'undefined') {
  window.busuanzi = busuanzi;
}
