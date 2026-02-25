# Coding for Life Sciences

## 建站信息

- **初建时间**: 2018年8月
- **重建时间**: 2020年4月
- **建站工具**: [Hexo](https://hexo.io/)
- **主题**: 基于 [jsimple](https://github.com/SumiMakito/hexo-theme-jsimple) 主题，并进行了自定义修改。

---

## 维护指南

### 主题修改与同步
若要保留对主题的修改，请先删除主题目录下的 `.git` 文件夹，然后运行：
```shell
git rm --cached themes/jsimple
git add themes/jsimple
git commit -m "Save modified theme"
```

### 依赖包更新
#### 常规更新
```shell
npm update
npm audit fix
git add .
git commit -m "update depend package"
git push origin source
```

#### 交互式升级 (推荐)
使用 `npm-upgrade` 交互式选择并升级 `package.json` 中的依赖版本：
```shell
npm install -g npm-upgrade
npm-upgrade
```

### 常见问题 (Troubleshooting)
- **RMarkdown 支持 (2020-04-10)**: 使用 blogdown 时，code trunk 的图片和缓存路径可能需要手动调整（推测为原主题 Bug）。
- **`spawnSync pandoc ETIMEDOUT` 错误 (2026-02-25)**: 由于 `hexo-renderer-pandoc` 默认超时时间较短（5000ms），在渲染包含大量复杂 LaTeX 公式的 Markdown 文件时易触发超时中断。解决方案：在根目录 `_config.yml` 中显式配置 `pandoc: \n  timeout: 60000` 将超时时间延长至 60 秒。

---

## 使用方法

### 发布新博文
```shell
hexo n post "title"
hexo g -d
git add . 
git commit -m 'publish new'
git push origin source
```

### 本地调试
如需使用本地安装的 Hexo：
```shell
./node_modules/hexo/bin/hexo g
```

---

## 更新日志

- **2026-02-25**: 升级依赖包，修复 `minimatch`。引入 `npm-upgrade` 管理依赖。修复 `hexo-renderer-pandoc` 的 `spawnSync pandoc ETIMEDOUT` 超时构建问题。集成 `busuanzi.simple.js` 实现基于本地存储的 PV/UV 访问统计功能。
- **2025-01-12**: 增加本地 Hexo 调用方式说明。
- **2022-06-21**: 手动修复 `package.json` 和 `package-lock.json` 版本冲突。