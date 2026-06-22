import os
import subprocess
import itertools
from datetime import datetime

# === НАСТРОЙКИ ПЛАТФОРМ ===
# Укажите правильные пути к вашим исполняемым файлам
EXECUTABLES = {
    "GPU": "./build/md_simulation_cuda",
    "CPU": "./build/md_simulation_cpu"  # Убедитесь, что это имя совпадает с вашей CPU-сборкой
}

# Базовая папка для результатов этой серии запусков
BASE_OUT_DIR = f"experiments_dual_{datetime.now().strftime('%Y%m%d_%H%M%S')}"

# Сетки физических параметров
PARAM_GRID = {
    "-N": [4000, 10976],                # Число частиц
    "-T": [87.3, 90.0],                 # Температура (K)
    "-steps": [40000],                  # Количество шагов
    "-tequil": [10000],                 # Время релаксации
}

# Дополнительные флаги без значений (например, включение термостата)
EXTRA_FLAGS = ["-Noz"] 

# Список платформ, на которых нужно провести расчеты
PLATFORMS = ["GPU", "CPU"]

def run_simulation():
    os.makedirs(BASE_OUT_DIR, exist_ok=True)
    
    # Генерируем комбинации физических параметров
    keys, values = zip(*PARAM_GRID.items())
    phys_experiments = [dict(zip(keys, v)) for v in itertools.product(*values)]
    
    # Создаем полный список задач: каждая физическая комбинация * каждая платформа
    all_runs = list(itertools.product(phys_experiments, PLATFORMS))
    total_runs = len(all_runs)
    
    print(f"=== Найдено {len(phys_experiments)} физ. комбинаций.")
    print(f"=== Общее количество запусков (с учетом CPU/GPU): {total_runs} ===")
    
    for idx, (exp, platform) in enumerate(all_runs, 1):
        executable = EXECUTABLES[platform]
        
        # Проверяем, существует ли бинарник перед запуском
        if not os.path.exists(executable):
            print(f"[{idx}/{total_runs}] [ПРОПУСК] Исполняемый файл для {platform} не найден по пути {executable}")
            continue
            
        # Формируем имя папки для конкретного эксперимента (добавляем префикс платформы)
        param_str = "_".join([f"{k.strip('-')}{v}" for k, v in exp.items()])
        folder_name = f"run_{idx}_{platform}_{param_str}"
        run_dir = os.path.join(BASE_OUT_DIR, folder_name)
        os.makedirs(run_dir, exist_ok=True)
        
        # Строим аргументы командной строки
        args = [executable]
        for param, val in exp.items():
            args.extend([param, str(val)])
        args.extend(EXTRA_FLAGS)
        
        print(f"[{idx}/{total_runs}] Запуск [{platform}]: {' '.join(args)}")
        
        stdout_log_path = os.path.join(run_dir, "simulation.log")
        stderr_log_path = os.path.join(run_dir, "error.log")
        
        start_time = datetime.now()
        with open(stdout_log_path, "w") as fout, open(stderr_log_path, "w") as ferr:
            try:
                process = subprocess.run(
                    args, 
                    stdout=fout, 
                    stderr=ferr, 
                    text=True,
                    check=True
                )
                status = "SUCCESS"
            except subprocess.CalledProcessError as e:
                status = f"FAILED (code {e.returncode})"
                print(f"  -> [ОШИБКА] Симуляция завершилась неудачно. См. {stderr_log_path}")
            except Exception as e:
                status = f"CRITICAL ERROR ({str(e)})"
                print(f"  -> [КРИТИЧЕСКАЯ ОШИБКА]: {e}")
                
        end_time = datetime.now()
        duration = (end_time - start_time).total_seconds()
        
        # Перемещаем файл траектории, если он создался в корне
        local_trajectory = "trajectory.csv"
        if os.path.exists(local_trajectory):
            os.rename(local_trajectory, os.path.join(run_dir, f"trajectory_{platform}_{param_str}.csv"))

        # Запись в сводную таблицу (summary.csv)
        with open(os.path.join(BASE_OUT_DIR, "summary.csv"), "a") as summary:
            if summary.tell() == 0:
                header = "Platform," + ",".join(keys) + ",Status,Duration_Sec,Folder\n"
                summary.write(header)
            
            row = f"{platform}," + ",".join([str(exp[k]) for k in keys]) + f",{status},{duration:.2f},{folder_name}\n"
            summary.write(row)

    print(f"\n=== Все расчеты завершены! Результаты сохранены в: {BASE_OUT_DIR} ===")

if __name__ == "__main__":
    run_simulation()