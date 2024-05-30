# zmp-network
Correlation network of zebrafish development

```
scratch 
cd detct/grcz11/
export gitdir=~/checkouts/zmp-network
export basedir=$(pwd)
screen -S zmp-network
```

Copy data
```
mkdir -p everything/filter-strict/
rclone copy --progress --max-depth 1 --filter "+ all.csv.gz" --filter "- *" sharepoint-qmul-buschlab:detct/grcz11/everything/filter-strict/ everything/filter-strict/
```

