# Coding for Life Sciences

建站信息：

    初建于2018年8月，重建于2020年4月。

    建站工具：hexo

    基于jsimple主题，有修改。


+ 保留修改后的主题的方法：删除主题中的.git隐藏文件夹，使用git rm --cached themes\jsimple; git add themes\jsimple;然后git commit。

+ issue:20200410 添加rmarkdown的方法，blogdown, code trunk图片和缓存路径需要手动改。。。目测是原作者的bug...

+ issue:20211027 日常更新与升级hexo依赖的javascript包：

```
    npm update
    npm audit fix
    git add .
    git commit -m "update depend package"
    git push origin source
```
Publish new content
```
    hexo n post "title"
    hexo g -d
    git add . 
    git commit -m 'publish new'
    git push origin source
```

20220621
manually fix, change version of package.json file, package-lock.json